#include <iostream>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>

using namespace std;

#define max(A, B) ((A) < (B)) ? (A) : (B);
const int NUM_THREADS = 4;

//===== Globalus kintamieji ===================================================

const int numDP = 5000; // Vietoviu skaicius (demand points, max 10000)
const int numPF = 5;	// Esanciu objektu skaicius (preexisting facilities)
const int numCL = 50;	// Kandidatu naujiems objektams skaicius (candidate locations)
const int numX = 3;		// Nauju objektu skaicius

double **demandPoints;	 // Geografiniai duomenys
double **distanceMatrix; // Masyvas atstumu matricai saugoti

int *solution = new int[numX];		// Naujas sprendinys
int *best_solution = new int[numX]; // Geriausias rastas sprendinys
double utility, best_utility;		// Naujo sprendinio ir geriausio sprendinio naudingumas (utility)

//===== Funkciju prototipai ===================================================

double getTime();
void loadDemandPoints();
double HaversineDistance(double lat1, double lon1, double lat2, double lon2);
double HaversineDistance(int i, int j);
double evaluateSolution(int *solution);
int increaseX(int *solution, int index, int maxindex);
void display_results(char *filename);

//=============================================================================

int main() {
    loadDemandPoints();			// Duomenu nuskaitymas is failo
    double t_start = getTime(); // Algoritmo vykdymo pradzios laikas

    //----- Atstumu matricos skaiciavimas -------------------------------------
    distanceMatrix = new double *[numDP];
    for (int i = 0; i < numDP; i++)
    {
        distanceMatrix[i] = new double[i + 1];

        for (int j = 0; j <= i; j++)
        {
            distanceMatrix[i][j] = HaversineDistance(
                demandPoints[i][0],
                demandPoints[i][1],
                demandPoints[j][0],
                demandPoints[j][1]);
        }
    }
    double t_matrix = getTime();
    cout << "Matricos skaiciavimo trukme: " << t_matrix - t_start << endl;

    //----- Pradines naujo ir geriausio sprendiniu reiksmes -------------------
    for (int i = 0; i < numX; i++)
    { // Pradines naujo ir geriausio sprendiniu koordinates: [0,1,2,...]
        solution[i] = i;
        best_solution[i] = i;
    }
    utility = evaluateSolution(solution); // Naujo sprendinio naudingumas (utility)
    best_utility = utility;				  // Geriausio sprendinio sprendinio naudingumas

    double t_new = getTime();
    cout << "Naujo sprendinio skaiciavimo trukme: " << t_new - t_matrix << endl;

    //----- Visų galimų sprendinių perrinkimas --------------------------------
    omp_set_num_threads(NUM_THREADS);
    bool increased = true;

    while (increased) {
        increased = increaseX(solution, numX - 1, numCL);
    
        utility = evaluateSolution(solution);
    
        if (utility > best_utility) {
            best_utility = utility;
            memcpy(best_solution, solution, sizeof(int) * numX);
        }
    }

    //----- Rezultatu spausdinimas --------------------------------------------
    double t_finish = getTime(); // Skaiciavimu pabaigos laikas
    cout << "Sprendinio paieskos trukme: " << t_finish - t_matrix << endl;
    cout << "Algoritmo vykdymo trukme: " << t_finish - t_start << endl;

    display_results("stdout");
    display_results("new.dat");
}

//===== Funkciju implementacijos (siu funkciju LYGIAGRETINTI NEREIKIA) ========

void loadDemandPoints()
{
    FILE *f;
    f = fopen("demandPoints.dat", "r");
    demandPoints = new double *[numDP];
    for (int i = 0; i < numDP; i++)
    {
        demandPoints[i] = new double[3];
        fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
    }
    fclose(f);
}

//=============================================================================

double HaversineDistance(double lat1, double lon1, double lat2, double lon2)
{
    double dlat = fabs(lat1 - lat2);
    double dlon = fabs(lon1 - lon2);
    double aa = pow((sin((double)dlat / (double)2 * 0.01745)), 2) + cos(lat1 * 0.01745)
                * cos(lat2 * 0.01745) * pow((sin((double)dlon / (double)2 * 0.01745)), 2);
    double c = 2 * atan2(sqrt(aa), sqrt(1 - aa));
    double d = 6371 * c;
    return d;
}

double HaversineDistance(int i, int j)
{
    if (i >= j)
        return distanceMatrix[i][j];
    else
        return distanceMatrix[j][i];
}

//=============================================================================

double getTime()
{
    struct timeval laikas;
    gettimeofday(&laikas, NULL);
    double rez = (double)laikas.tv_sec + (double)laikas.tv_usec / 1000000;
    return rez;
}

//=============================================================================

/* Mutates nothing */
double evaluateSolution(int *solution)
{
    double _utility = 0;
    double _total_utility = 0;
    int _bestPF;
    int _best_solution;
    double _d;

    #pragma omp parallel reduction (+:_utility, _total_utility) private(_bestPF, _best_solution, _d)
    #pragma omp for
    for (int i = 0; i < numDP; i++) {
        _bestPF = 1e5;

        for (int j = 0; j < numPF; j++)
        {
            _d = HaversineDistance(i, j);
            _bestPF = max(_d, _bestPF);
        }

        _best_solution = 1e5;
        for (int j = 0; j < numX; j++)
        {
            _d = HaversineDistance(i, solution[j]);
            _best_solution = max(_d, _best_solution);
        }

        _utility += (_best_solution < _bestPF)  ?       demandPoints[i][2] : 0;
        _utility += (_best_solution == _bestPF) ? 0.3 * demandPoints[i][2] : 0;

        _total_utility += demandPoints[i][2];
    }
    return _utility / _total_utility * 100;
}

//=============================================================================

/* Mutates: solution */      /* numX=3 - 1    numCL=50 */
int increaseX(int *solution, int index, int maxindex)
{
    if (solution[index] + 1 < maxindex - (numX - index - 1))
    {
        solution[index]++;
    }
    else
    {
        if ((index == 0) && (solution[index] + 1 == maxindex - (numX - index - 1)))
        {
            return 0;
        }
        else
        {
            if (increaseX(solution, index - 1, maxindex))
            {
                solution[index] = solution[index - 1] + 1;
            }
            else
                return 0;
        }
    }
    return 1;
}

void display_results(char *filename)
{
    const char *cmp = "stdout";
    FILE *fp;

    if (strcmp(filename, cmp) == 0) {
        fp = stdout;
    } else {
        fp = fopen(filename, "w+");
    }

    for (int i = 0; i < numX; i++) { fprintf(fp, "%i\t", best_solution[i]); }
    fprintf(fp, "\t%.3f\n", best_utility);
}
