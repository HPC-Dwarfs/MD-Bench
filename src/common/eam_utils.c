/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <allocate.h>
#include <atom.h>
#include <eam.h>
#include <parameter.h>
#include <util.h>

#ifndef MAXLINE
#define MAXLINE 4096
#endif
Eam eam;

void initEam(Parameter* param)
{
    int ntypes = param->ntypes;
    eam.nmax   = 0;
    eam.fp     = NULL;

    readEamFile(&eam.file, param->eam_file);
    param->mass     = eam.file.mass;
    param->cutforce = eam.file.cut;
    param->cutneigh = param->cutforce + 1.0;
    param->temp     = 600.0;
    param->dt       = 0.001;
    param->rho      = 0.07041125;
    param->dtforce  = 0.5 * param->dt / param->mass;

    // convert read-in file(s) to arrays and spline them
    file2array(&eam);
    array2spline(&eam, param);
}

void readEamFile(Funcfl* file, const char* filename)
{
    FILE* fptr;
    char line[MAXLINE];

    fptr = fopen(filename, "r");
    if (fptr == NULL) {
        printf("Can't open EAM Potential file: %s\n", filename);
        exit(0);
    }

    int tmp;
    readline(line, fptr);
    readline(line, fptr);
#if PRECISION == 2
    sscanf(line, "%d %lg", &tmp, &(file->mass));
    readline(line, fptr);
    sscanf(line,
        "%d %lg %d %lg %lg",
        &file->nrho,
        &file->drho,
        &file->nr,
        &file->dr,
        &file->cut);
#else
    sscanf(line, "%d %g", &tmp, &(file->mass));
    readline(line, fptr);
    sscanf(line,
        "%d %g %d %g %g",
        &file->nrho,
        &file->drho,
        &file->nr,
        &file->dr,
        &file->cut);
#endif

    // printf("Read: %lf %i %lf %i %lf
    // %lf\n",file->mass,file->nrho,file->drho,file->nr,file->dr,file->cut);
    file->frho = (MD_FLOAT*)allocate(ALIGNMENT, (file->nrho + 1) * sizeof(MD_FLOAT));
    file->rhor = (MD_FLOAT*)allocate(ALIGNMENT, (file->nr + 1) * sizeof(MD_FLOAT));
    file->zr   = (MD_FLOAT*)allocate(ALIGNMENT, (file->nr + 1) * sizeof(MD_FLOAT));
    grab(fptr, file->nrho, file->frho);
    grab(fptr, file->nr, file->zr);
    grab(fptr, file->nr, file->rhor);
    for (int i = file->nrho; i > 0; i--)
        file->frho[i] = file->frho[i - 1];
    for (int i = file->nr; i > 0; i--)
        file->rhor[i] = file->rhor[i - 1];
    for (int i = file->nr; i > 0; i--)
        file->zr[i] = file->zr[i - 1];
    fclose(fptr);
}

void file2array(Eam* eam)
{
    int i, j, k, m, n;
    double sixth = 1.0 / 6.0;

    // determine max function params from all active funcfl files
    // active means some element is pointing at it via map
    int active;
    double rmax, rhomax;
    eam->dr = eam->drho = rmax = rhomax = 0.0;
    active                              = 0;
    Funcfl* file                        = &eam->file;
    eam->dr                             = MAX(eam->dr, file->dr);
    eam->drho                           = MAX(eam->drho, file->drho);
    rmax                                = MAX(rmax, (file->nr - 1) * file->dr);
    rhomax                              = MAX(rhomax, (file->nrho - 1) * file->drho);

    // set nr,nrho from cutoff and spacings
    // 0.5 is for round-off in divide
    eam->nr   = (int)(rmax / eam->dr + 0.5);
    eam->nrho = (int)(rhomax / eam->drho + 0.5);

    // ------------------------------------------------------------------
    // setup frho arrays
    // ------------------------------------------------------------------

    // allocate frho arrays
    // nfrho = # of funcfl files + 1 for zero array
    eam->frho = (MD_FLOAT*)allocate(ALIGNMENT, (eam->nrho + 1) * sizeof(MD_FLOAT));

    // interpolate each file's frho to a single grid and cutoff
    double r, p, cof1, cof2, cof3, cof4;
    n = 0;

    for (m = 1; m <= eam->nrho; m++) {
        r = (m - 1) * eam->drho;
        p = r / file->drho + 1.0;
        k = (int)(p);
        k = MIN(k, file->nrho - 2);
        k = MAX(k, 2);
        p -= k;
        p            = MIN(p, 2.0);
        cof1         = -sixth * p * (p - 1.0) * (p - 2.0);
        cof2         = 0.5 * (p * p - 1.0) * (p - 2.0);
        cof3         = -0.5 * p * (p + 1.0) * (p - 2.0);
        cof4         = sixth * p * (p * p - 1.0);
        eam->frho[m] = cof1 * file->frho[k - 1] + cof2 * file->frho[k] +
                       cof3 * file->frho[k + 1] + cof4 * file->frho[k + 2];
    }

    // ------------------------------------------------------------------
    // setup rhor arrays
    // ------------------------------------------------------------------

    // allocate rhor arrays
    // nrhor = # of funcfl files
    eam->rhor = (MD_FLOAT*)allocate(ALIGNMENT, (eam->nr + 1) * sizeof(MD_FLOAT));

    // interpolate each file's rhor to a single grid and cutoff
    for (m = 1; m <= eam->nr; m++) {
        r = (m - 1) * eam->dr;
        p = r / file->dr + 1.0;
        k = (int)(p);
        k = MIN(k, file->nr - 2);
        k = MAX(k, 2);
        p -= k;
        p            = MIN(p, 2.0);
        cof1         = -sixth * p * (p - 1.0) * (p - 2.0);
        cof2         = 0.5 * (p * p - 1.0) * (p - 2.0);
        cof3         = -0.5 * p * (p + 1.0) * (p - 2.0);
        cof4         = sixth * p * (p * p - 1.0);
        eam->rhor[m] = cof1 * file->rhor[k - 1] + cof2 * file->rhor[k] +
                       cof3 * file->rhor[k + 1] + cof4 * file->rhor[k + 2];
        // if(m==119)printf("BuildRho: %e %e %e %e %e
        // %e\n",rhor[m],cof1,cof2,cof3,cof4,file->rhor[k]);
    }

    // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
    // for funcfl files, I,J mapping only depends on I
    // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

    // ------------------------------------------------------------------
    // setup z2r arrays
    // ------------------------------------------------------------------

    // allocate z2r arrays
    // nz2r = N*(N+1)/2 where N = # of funcfl files
    eam->z2r = (MD_FLOAT*)allocate(ALIGNMENT, (eam->nr + 1) * sizeof(MD_FLOAT));

    // create a z2r array for each file against other files, only for I >= J
    // interpolate zri and zrj to a single grid and cutoff
    double zri, zrj;
    Funcfl* ifile = &eam->file;
    Funcfl* jfile = &eam->file;

    for (m = 1; m <= eam->nr; m++) {
        r = (m - 1) * eam->dr;
        p = r / ifile->dr + 1.0;
        k = (int)(p);
        k = MIN(k, ifile->nr - 2);
        k = MAX(k, 2);
        p -= k;
        p    = MIN(p, 2.0);
        cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
        cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
        cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
        cof4 = sixth * p * (p * p - 1.0);
        zri  = cof1 * ifile->zr[k - 1] + cof2 * ifile->zr[k] + cof3 * ifile->zr[k + 1] +
              cof4 * ifile->zr[k + 2];

        p = r / jfile->dr + 1.0;
        k = (int)(p);
        k = MIN(k, jfile->nr - 2);
        k = MAX(k, 2);
        p -= k;
        p    = MIN(p, 2.0);
        cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
        cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
        cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
        cof4 = sixth * p * (p * p - 1.0);
        zrj  = cof1 * jfile->zr[k - 1] + cof2 * jfile->zr[k] + cof3 * jfile->zr[k + 1] +
              cof4 * jfile->zr[k + 2];

        eam->z2r[m] = 27.2 * 0.529 * zri * zrj;
    }
}

void array2spline(Eam* eam, Parameter* param)
{
    eam->rdr      = 1.0 / eam->dr;
    eam->rdrho    = 1.0 / eam->drho;
    eam->nrho_tot = (eam->nrho + 1) * 7 + 64;
    eam->nr_tot   = (eam->nr + 1) * 7 + 64;
    eam->nrho_tot -= eam->nrho_tot % 64;
    eam->nr_tot -= eam->nr_tot % 64;

    int ntypes       = param->ntypes;
    eam->frho_spline = (MD_FLOAT*)allocate(ALIGNMENT,
        ntypes * ntypes * eam->nrho_tot * sizeof(MD_FLOAT));
    eam->rhor_spline = (MD_FLOAT*)allocate(ALIGNMENT,
        ntypes * ntypes * eam->nr_tot * sizeof(MD_FLOAT));
    eam->z2r_spline  = (MD_FLOAT*)allocate(ALIGNMENT,
        ntypes * ntypes * eam->nr_tot * sizeof(MD_FLOAT));
    interpolate(eam->nrho, eam->drho, eam->frho, eam->frho_spline);
    interpolate(eam->nr, eam->dr, eam->rhor, eam->rhor_spline);
    interpolate(eam->nr, eam->dr, eam->z2r, eam->z2r_spline);

    // replicate data for multiple types;
    for (int tt = 0; tt < ntypes * ntypes; tt++) {
        for (int k = 0; k < eam->nrho_tot; k++)
            eam->frho_spline[tt * eam->nrho_tot + k] = eam->frho_spline[k];
        for (int k = 0; k < eam->nr_tot; k++)
            eam->rhor_spline[tt * eam->nr_tot + k] = eam->rhor_spline[k];
        for (int k = 0; k < eam->nr_tot; k++)
            eam->z2r_spline[tt * eam->nr_tot + k] = eam->z2r_spline[k];
    }
}

void interpolate(int n, MD_FLOAT delta, MD_FLOAT* f, MD_FLOAT* spline)
{
    for (int m = 1; m <= n; m++)
        spline[m * 7 + 6] = f[m];

    spline[1 * 7 + 5]       = spline[2 * 7 + 6] - spline[1 * 7 + 6];
    spline[2 * 7 + 5]       = 0.5 * (spline[3 * 7 + 6] - spline[1 * 7 + 6]);
    spline[(n - 1) * 7 + 5] = 0.5 * (spline[n * 7 + 6] - spline[(n - 2) * 7 + 6]);
    spline[n * 7 + 5]       = spline[n * 7 + 6] - spline[(n - 1) * 7 + 6];

    for (int m = 3; m <= n - 2; m++)
        spline[m * 7 + 5] = ((spline[(m - 2) * 7 + 6] - spline[(m + 2) * 7 + 6]) +
                                8.0 *
                                    (spline[(m + 1) * 7 + 6] - spline[(m - 1) * 7 + 6])) /
                            12.0;

    for (int m = 1; m <= n - 1; m++) {
        spline[m * 7 + 4] = 3.0 * (spline[(m + 1) * 7 + 6] - spline[m * 7 + 6]) -
                            2.0 * spline[m * 7 + 5] - spline[(m + 1) * 7 + 5];
        spline[m * 7 + 3] = spline[m * 7 + 5] + spline[(m + 1) * 7 + 5] -
                            2.0 * (spline[(m + 1) * 7 + 6] - spline[m * 7 + 6]);
    }

    spline[n * 7 + 4] = 0.0;
    spline[n * 7 + 3] = 0.0;

    for (int m = 1; m <= n; m++) {
        spline[m * 7 + 2] = spline[m * 7 + 5] / delta;
        spline[m * 7 + 1] = 2.0 * spline[m * 7 + 4] / delta;
        spline[m * 7 + 0] = 3.0 * spline[m * 7 + 3] / delta;
    }
}

void grab(FILE* fptr, int n, MD_FLOAT* list)
{
    char* ptr;
    char line[MAXLINE];
    int i = 0;

    while (i < n) {
        readline(line, fptr);
        ptr       = strtok(line, " \t\n\r\f");
        list[i++] = atof(ptr);
        while ((ptr = strtok(NULL, " \t\n\r\f")))
            list[i++] = atof(ptr);
    }
}
