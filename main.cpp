#include <iostream>
#include <iomanip>
#include <cmath> 
#include <cstdio>   
#include <sys/stat.h> 
#include <sys/types.h> 
#include <string>   
#include <cstdlib>

using namespace std;


const int N = 51;
const double dx = 1.0 / (N-1);
const double dy = 1.0 / (N-1);
const double cfl = 0.5;
const double dt = dx * cfl;//0.005;
const int iterations = 4000;
const double Re = 1000;
const int output_interval = 10;
const int possion_iterations = 10000;
double poisson_tolerance = 1e-6;


void initialize(double s[N][N], double w[N][N], double u[N][N],double v[N][N]);             // initialize fields
void computeVorticity(double s[N][N], double w[N][N], double u[N][N],double v[N][N]);       // compute vorticity
void computeStreamFunction(double s[N][N], double w[N][N]);                                 // compute stream function
void computeVelocity(double s[N][N], double u[N][N],double v[N][N]);                        // compute velocities
void display(double data[N][N]);                                                            // display matrix
void saveVTK(const std::string& fieldName, const double data[N][N]);                        // store results


int main() {
    double s[N][N];
    double w[N][N];
    double u[N][N];
    double v[N][N];
    initialize(s,w,u,v);
    for (int t=0; t<iterations; t++){
        printf("iteration:  %d     ",t);
        computeVorticity(s,w,u,v);

        computeStreamFunction(s,w);

        computeVelocity(s,u,v);

        if (t%output_interval == 0){
            saveVTK("vorticity", w);
            saveVTK("streamFunc", s);
            saveVTK("u", u);
            saveVTK("v", v);
            //display(u);
        }
    }

    return 0;
}




// function to initialize fields
void initialize(double s[N][N], double w[N][N], double u[N][N],double v[N][N]){
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            s[i][j] = 0;
            w[i][j] = 0;
            u[i][j] = 0;
            v[i][j] = 0;
            if (i==N-1){
                u[i][j] = 1;
            }
        }
    }  

}

// function to compute Vorticity
void computeVorticity(double s[N][N], double w[N][N], double u[N][N],double v[N][N]){
    // apply BC for w
    for (int i=0; i<N; i++){
        // bottom and top
        w[0][i] = -2 * s[1][i]/(dy*dy);
        w[N-1][i] = -2 * (s[N-2][i]/(dy*dy) + 1/dy);
        // left and right
        w[i][0] = -2 * s[i][1]/(dx*dx);
        w[i][N-1] = -2 * s[i][N-2]/(dx*dx);
        
    }
    
    // copy w to w0
    double w0[N][N]; 
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            w0[i][j] = w[i][j];
        }
    }
     
    // transport eq for w
    for (int i=1; i<N-1; i++){
        for (int j=1; j<N-1; j++){
            w[i][j] = w0[i][j] - dt * (
                                        (u[i][j] * ((w0[i][j+1] - w0[i][j-1]) / (2*dx))) +
                                        (v[i][j] * ((w0[i+1][j] - w0[i-1][j]) / (2*dy))) -
                                        ((1./Re) * ((w0[i][j+1]- 2*w0[i][j]+ w0[i][j-1]) / (dx*dx))) -
                                        ((1./Re) * ((w0[i+1][j]- 2*w0[i][j]+ w0[i-1][j]) / (dy*dy)))
                                      );
        }
    }
    
}

// function to compute stream function
void computeStreamFunction(double s[N][N], double w[N][N]){
    
    double s0[N][N];
    for(int t=0; t<possion_iterations; t++){
        //copy s to s0
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                s0[i][j] = s[i][j];
            }
        }
        // poisson eq for s
        for (int i=1; i<N-1; i++){
            for (int j=1; j<N-1; j++){
                s[i][j] = (
                            (s0[i][j+1] + s0[i][j-1]) * (dy*dy) + 
                            (s0[i+1][j] + s0[i-1][j]) * (dx*dx) +
                            dx*dx*dy*dy*w[i][j]
                          ) /  (2*(dx*dx + dy*dy));
            }
        }
        double e = 0.0;
        for (int i=1; i < N-1; i++)
        {
            for (int j=1; j < N-1; j++)
            {
            e += fabs(s[i][j] - s0[i][j]);
            }
        }
        //e /= (N-1)*(N-1);
        if(e<poisson_tolerance){
            printf("poisson equation iteration: %d \n", t);
            return;
        }
    }
    printf("Error: max number of iteration for posson is achieved");
    exit(1);
     
}

// function to compute velocities
void computeVelocity(double s[N][N], double u[N][N],double v[N][N]){
    for (int i=1; i<N-1; i++){
            for (int j=1; j<N-1; j++){
                u[i][j] =  (s[i+1][j] - s[i-1][j])/(2*dy);
                v[i][j] = -(s[i][j+1] - s[i][j-1])/(2*dx);
            }
    }
}

// function to display matric
void display(double data[N][N]){
    cout << fixed << setprecision(4);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            cout << setw(8) << data[i][j];
        }
        cout << endl; 
    }
    cout << endl;
}

// function to store results
void saveVTK(const std::string& fieldname, const double data[N][N]) {
    char name[64];
    FILE *pf;
    static int count = 0;

    // Ensure output directory exists
    struct stat st = {0};
    if (stat("./results_vtk", &st) == -1) {
        mkdir("./results_vtk", 0700);
    }

    // Create filename using count
    snprintf(name, sizeof(name), "./results_vtk/%s-%03d.vtk", fieldname.c_str(), count);

    // Open the file for writing
    pf = fopen(name, "w");
    if (pf == nullptr) {
        perror("Error while opening file");
        return;
    }

    // Write VTK file header
    fprintf(pf, "# vtk DataFile Version 2.0\n");
    fprintf(pf, "Field data for %s\n", fieldname.c_str());
    fprintf(pf, "ASCII\n");
    fprintf(pf, "DATASET STRUCTURED_POINTS\n");
    fprintf(pf, "DIMENSIONS %d %d 1\n", N, N);
    fprintf(pf, "ORIGIN 0 0 0\n");
    fprintf(pf, "SPACING 1 1 1\n");
    fprintf(pf, "POINT_DATA %d\n", N * N);
    fprintf(pf, "SCALARS %s float 1\n", fieldname.c_str());
    fprintf(pf, "LOOKUP_TABLE default\n");

    // Write the field data
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(pf, "%.6f ", data[i][j]);
        }
        fprintf(pf, "\n");
    }

    // Close the file
    fclose(pf);
    count++;  // Increment the file count for the next call
}
