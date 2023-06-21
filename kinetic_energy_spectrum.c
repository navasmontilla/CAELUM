/*
    Para compilar: gcc cascade.c -o cascade.exe -I kissfft kissfft/kiss_fftnd.c -I kissfft kissfft/kiss_fft.c

    Usa la libreria kissfft (1.3.0): https://github.com/mborgerding/kissfft
*/ 
    

#include "kissfft/kiss_fftnd.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define fran (rand()/((double)RAND_MAX+1))
#define frand (fran*2-1)
#define ABS(x) (x < 0 ? -x : x)


double* fft2d_magnitude(double* img, double* mag, int N);
double* fft3d_magnitude(double* img, double* mag, int N);
void compute_energy_cascade_2D(double* U, double* V, double* EK, int N, int box_radius);
void compute_energy_cascade_3D(double* U, double* V, double* W, double* EK, int N, int box_radius);
void copy_array(double* a, double* b, int N);
void fftshift_2D(double* orig, int N);
void fftshift_3D(double* orig, int N);

int write_column(char* fname, double *data, int N);
int read_file_list_2D(char* fname, double* U, double* V, int N);
int read_file_list_3D(char* fname, double* U, double* V, double* W, int N);
void read_file_column(char* fname, double* arr, int N);


double* fft2d_magnitude(double* img, double* mag, int N){

    int dims[] = {N, N};
    kiss_fftnd_cfg config = kiss_fftnd_alloc(dims, 2, 0, NULL, NULL);

    kiss_fft_cpx* in  = malloc(N * N * sizeof(kiss_fft_cpx));
    kiss_fft_cpx* out = malloc(N * N * sizeof(kiss_fft_cpx));    

    // Llenamos los datos de entrada con img
    for(int i = 0; i < N*N; i++){
        in[i].r = img[i];
        in[i].i = 0.0;
    }

    // realiza la transformada de Fourier 2D
    kiss_fftnd(config, (kiss_fft_cpx*)in, (kiss_fft_cpx*)out);

    // calcula la magnitud de la transformada
    for(int i = 0; i < N*N; i++){
        mag[i] = sqrt(out[i].r*out[i].r + out[i].i*out[i].i);
        mag[i] = mag[i] / (N*N);
    }

    // libera la memoria
    // free(config);
    // free(in);
    free(out);
    return mag;
}

double* fft3d_magnitude(double* img, double* mag, int N){

    int i;
    int dims[] = {N, N, N};
    kiss_fftnd_cfg config = kiss_fftnd_alloc(dims, 3, 0, NULL, NULL);

    kiss_fft_cpx* in  = malloc(N * N * N * sizeof(kiss_fft_cpx));
    kiss_fft_cpx* out = malloc(N * N * N * sizeof(kiss_fft_cpx));

    // Llenamos los datos de entrada con img
    for(int i = 0; i < N*N*N; i++){
        in[i].r = img[i];
        in[i].i = 0.0;
    }

    // realiza la transformada de Fourier 2D
    kiss_fftnd(config, (kiss_fft_cpx*)in, (kiss_fft_cpx*)out);

    // calcula la magnitud de la transformada
    for(i = 0; i < N*N*N; i++){
        mag[i] = sqrt(out[i].r*out[i].r + out[i].i*out[i].i) / (N*N*N);
    }

    // libera la memoria
    // free(config);
    // free(in);
    free(out);
    return mag;
}

void fftshift_2D(double* orig, int N){
    // Calcular los desplazamientos
    int i_shift, j_shift, i, j;
    int shift = N / 2;

    double* shifted = malloc(N * N * sizeof(double));
    
    // Copiar los elementos de la matriz de entrada a la matriz de salida con el desplazamiento adecuado
    for(int i = 0; i < N; i++){
        i_shift = (i + shift) % N;
        for(int j = 0; j < N; j++){
            j_shift = (j + shift) % N;
            shifted[i*N + j] = orig[i_shift*N + j_shift];
        }
    }

    copy_array(shifted, orig, N*N);
    free(shifted);
}

void fftshift_3D(double* orig, int N){
    // Calcular los desplazamientos
    int i_shift, j_shift, k_shift, i, j, k;
    int shift = N / 2;

    double* shifted = malloc(N * N * N * sizeof(double));
    
    // Copiar los elementos del arreglo de entrada al arreglo de salida con el desplazamiento adecuado
    for(int i = 0; i < N; i++){
        i_shift = (i + shift) % N;
        for(int j = 0; j < N; j++){
            j_shift = (j + shift) % N;
            for(int k = 0; k < N; k++){
                k_shift = (k + shift) % N;
                shifted[i*N*N + j*N + k] = orig[i_shift*N*N + j_shift*N + k_shift];
            }
        }
    }

    copy_array(shifted, orig, N*N*N);
    free(shifted);
}

void compute_energy_cascade_2D(double* U, double* V, double* EK, int N, int box_radius){

    int center = N/2;
    int i, j, dist;

    double* EK_U = malloc(N * N * sizeof(double));
    double* EK_V = malloc(N * N * sizeof(double));
    double* EK_U_avg = malloc(box_radius  * sizeof(double));
    double* EK_V_avg = malloc(box_radius  * sizeof(double));

    fft2d_magnitude(U, EK_U, N);
    fft2d_magnitude(V, EK_V, N);
    
    fftshift_2D(EK_U, N);
    fftshift_2D(EK_V, N);

    for (i = 0; i < N*N; i++){
        EK_U[i] = EK_U[i] * EK_U[i];
        EK_V[i] = EK_V[i] * EK_V[i];
    }

    for (i = 0; i < box_radius; i++){
        EK_U_avg[i] = 1E-50;
        EK_V_avg[i] = 1E-50;
    }

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            dist = sqrt((i-center)*(i-center) + (j-center)*(j-center));
            if (dist<box_radius){
                EK_U_avg[dist] += EK_U[i*N+j];
                EK_V_avg[dist] += EK_V[i*N+j];
            }
        }
    }

    for (i = 0; i < N/2+1; i++){
        EK[i] = 0.5 * (EK_U_avg[i] + EK_V_avg[i]);
    }

    free(EK_U);
    free(EK_V);
    free(EK_U_avg);
    free(EK_V_avg);  
}

void compute_energy_cascade_3D(double* U, double* V, double* W, double* EK, int N, int box_radius){

    int center = N/2;
    int i, j, k, dist;

    double* EK_U = malloc(N * N * N * sizeof(double));
    double* EK_V = malloc(N * N * N * sizeof(double));
    double* EK_W = malloc(N * N * N * sizeof(double));
    double* EK_U_sum = malloc(box_radius  * sizeof(double));
    double* EK_V_sum = malloc(box_radius  * sizeof(double));
    double* EK_W_sum = malloc(box_radius  * sizeof(double));

    fft3d_magnitude(U, EK_U, N);
    fft3d_magnitude(V, EK_V, N);
    fft3d_magnitude(W, EK_W, N);
    
    fftshift_3D(EK_U, N);
    fftshift_3D(EK_V, N);
    fftshift_3D(EK_W, N);

    for (i = 0; i < N*N*N; i++){
        EK_U[i] = EK_U[i] * EK_U[i];
        EK_V[i] = EK_V[i] * EK_V[i];
        EK_W[i] = EK_W[i] * EK_W[i];
    }

    for (i = 0; i < box_radius; i++){
        EK_U_sum[i] = 1E-50;
        EK_V_sum[i] = 1E-50;
        EK_W_sum[i] = 1E-50;
    }

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            for (k = 0; k < N; k++){
                dist = sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
                if (dist<box_radius){
                    EK_U_sum[dist] += EK_U[i*N*N + j*N + k];
                    EK_V_sum[dist] += EK_V[i*N*N + j*N + k];
                    EK_W_sum[dist] += EK_W[i*N*N + j*N + k];
                }
            }
        }
    }

    for (i = 0; i < N/2+1; i++){
        EK[i] = 0.5 * (EK_U_sum[i] + EK_V_sum[i] + EK_W_sum[i]);
    }

    free(EK_U);
    free(EK_V);
    free(EK_W);
    free(EK_U_sum);
    free(EK_V_sum);  
    free(EK_W_sum);  
}

int write_column(char* fname, double *data, int N){
    FILE* fp = fopen(fname, "w");

    if (fp == NULL){
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    for(int i = 0; i < N; i++){
        fprintf(fp, "%.14e\n", data[i]);
    }

    fclose(fp);

    return 0;
}

int read_file_list_2D(char* fname, double* U, double* V, int N){
    
    // Carga las velocidades U y V de un archivo de salida

    int i,j,idx;
    char line[110];
    double asd;

    FILE* fp = fopen(fname, "r");

    if (fp == NULL){
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    fscanf(fp, "%*[^\n]%*c");
    fscanf(fp, "%*[^\n]%*c");

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            idx = i*N + j;
            fscanf(fp, "%*f %*f %*f %le %le %*[^\n]%*c", &U[idx], &V[idx]);
            // printf("%le  %le\n", U[idx], V[idx]);
        }
    }

    fclose(fp);
    return 0;
}

int read_file_list_3D(char* fname, double* U, double* V, double* W, int N){

    // Carga las velocidades U, V y W de un archivo de salida

    int i,j,idx;
    char line[110];
    double asd;

    FILE* fp = fopen(fname, "r");

    if (fp == NULL){
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    fscanf(fp, "%*[^\n]%*c");
    fscanf(fp, "%*[^\n]%*c");

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            idx = i*N + j;
            fscanf(fp, "%*f %*f %*f %le %le %le%*[^\n]%*c", &U[idx], &V[idx], &W[idx]);
            // printf("%le  %le\n", U[idx], V[idx]);
        }
    }

    fclose(fp);
    return 0;
}

void read_file_column(char* fname, double* arr, int N){
    FILE* f = fopen(fname, "r");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            fscanf(f, "%le", &arr[i*N+j]);
        }
    }
    fclose(f);
}

void copy_array(double* a, double* b, int N){
    for (int i = 0; i < N; i++){
        b[i] = a[i];
    }
}

int main(){ 

    int Nc = 128;
    int boxR = Nc/2;
    char fname[100];

    printf("\n  N=%d    boxR = %d\n", Nc, boxR);

    printf("  Allocating memory...  ");
    double* U = malloc(Nc * Nc *  Nc * sizeof(double));
    double* V = malloc(Nc * Nc *  Nc * sizeof(double));
    double* W = malloc(Nc * Nc *  Nc * sizeof(double));
    double* EK = malloc(boxR * sizeof(double));
    printf("Done.\n");

    printf("  Reading files...  ");
    read_file_list_3D("list006.out", U, V, W, Nc);
    printf("Done.\n");

    printf("  Computing energy cascade...  ");
    compute_energy_cascade_3D(U, V, W, EK, Nc, boxR);
    printf("Done.\n");

    printf("  Printing file...  ");
    write_column("energy_spectrum_C.txt", EK, boxR-1);
    printf("Done.\n\n");

    free(U);
    free(V);
    free(W);
    free(EK);
    
    return 0;
}