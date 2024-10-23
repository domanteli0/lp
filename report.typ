#import "@preview/mitex:0.2.4": *
#import "@preview/plotst:0.2.0": *
#import "@preview/oxifmt:0.2.0": strfmt

#let pageref(label) = context {
  let loc = locate(label)
  let nums = counter(page).at(loc)
  link(label, str(nums.first()))
}

#set figure(
  supplement: [Pav.]
)

#show raw.where(block: true): it => box(
  fill: rgb("#e6e5e5"),
  outset: 3pt,
  radius: 2pt,
  width: 100%,
  align(left)[#it]
)

#align(center, text(18pt)[
  = *Laboratorinis darbai*
])

#align(center, text(15pt)[
  Domantas Keturakis \
  #text(13pt)[
    Spalis 2024
  ]
  ]
 )

#show heading.where(depth: 1): it => [
  #set text(17pt)
  #block(underline(smallcaps(it.body)))
]

#show heading.where(depth: 2): it => [
  #set text(15pt)
  #it.body
]

#set par(justify: true)

= Užduotis \#1

== Pirma dalis

=== Paralelizavimo galimybių analizė

Pirmasis skaičiavimo skaičiavimo ciklas (@first_cycle) palyginus užtrunka nedaug laiko (apie ~0.002 sekundės). Praktiniems tikslams, jį galima ignoruoti.

#figure(
  placement: none,
  ```c
  for (int i=0; i<numX; i++) {
      X[i] = i;
      bestX[i] = i;
  }
  u = evaluateSolution(X);
  bestU = u;
  ```,
  caption: [Pradinės naujo ir geriausio sprendinių reikšmių apskaičiavimas]
) <first_cycle>

Didžiąją dalį laiko užima šis ciklas:

#figure(
  placement: none,
```c
while (increaseX(X, numX-1, numCL) == true) {
    u = evaluateSolution(X);
    if (u > bestU) {
        bestU = u;
        for (int i=0; i<numX; i++) bestX[i] = X[i];
    }
}
```,
  caption: [Visų galimų sprendinių perrinkimas]
)

Šio ciklo iš esmės "aklai" parelilizuoti negalima, nes reikia atsižvelgti į tai kad:

  - `increaseX` - keičia masyvo `X` reikšmę (@increseX), ir ne tik `index`-ąjį elementą, rekursyviai kviesdamas save sumažina `index` reikšmę vienu, t.y. iškvietus `increaseX` visos masyvos reikšmės yra keičiamos. To pasekmė, kad `index`-ojo elemento skaičiavimo negalima paskirstyti skirtingoms gijoms, kitaip vėlesnėms gijoms reikėtų laukti, kol praeita gija baigs savo darbą, visiškai nustelbiant parelelizavimo naudą.

  - `increaseX` skaičiavimai priklauso vienas nuo kito, t.y. norint apskaičiuoti `X` reikšmę $n$-ame ciklo vykdyme, reikia pirma apskaičiuoti `X` reikšmę $(n-1)$-ame ciklo vykdyme. Analogiškai negalima paralelizuoti nes kitos gijos lauktų, kol praeita gija baigs savo darbą.

  - `if (u > bestU) { ... }` irgi gali tik vienas ciklas vienu metu, nes `bestU` ir `X` pakeitimas turi būti atliekamas "žingsniu" - t.y. atomiškai.

Iš esmės neperrašius `increaseX`, šios funkcijos ir jos kvietimo cikle, yra nepraktiška parelilizuoti.

#figure(
  placement: none,
```c
int increaseX(int *X, int index, int maxindex) {
    if (X[index]+1 < maxindex-(numX-index-1)) {
        X[index]++;
    }
    else {		 
        if ((index == 0) && (X[index]+1 == maxindex-(numX-index-1))) {
            return 0;
        }
        else {
            if (increaseX(X, index-1, maxindex)) X[index] = X[index-1]+1;
            else return 0;
        }	
    }
    return 1;
}
```,
  caption: [Funkcija `increseX`]
) <increseX>

Tuo tarpu funkcija `evaluateSolution` (@evaluateSolution) nekeičia jokių globalių kintamųjų ar savo argumentų. Analitiškai žiūrint galima spėti, kad čia ir didžioji dalis skaičiavimo laiko yra sugaištama. Teorinis šios funkcijos _big O_ yra $O(#raw("numDP") dot #raw("numX"))$, tuo tarpu `increseX` rekursyviai save gali iškviesti daugiausiai `numX` kartų.

#figure(
  placement: none,
```c
double evaluateSolution(int *X) {
    double U = 0;
    double totalU = 0;
    int bestPF;
    int bestX;
    double d;

    for (int i=0; i<numDP; i++) {
        totalU += demandPoints[i][2];
        bestPF = 1e5;
        for (int j=0; j<numPF; j++) {
            d = HaversineDistance(i, j);
            if (d < bestPF) bestPF = d;
        }
        bestX = 1e5;
        for (int j=0; j<numX; j++) {
            d = HaversineDistance(i, X[j]);
            if (d < bestX) bestX = d;
        }

        if (bestX < bestPF) U += demandPoints[i][2];
        else if (bestX == bestPF) U += 0.3*demandPoints[i][2];
    }
    return U/totalU*100;
}
```,
  caption: [Funkcija `evaluateSolution`]
) <evaluateSolution>


#pagebreak()
=== Sprendimo paralelizavimas

Dėl anksčiau išvardintų priežąsčių `increaseX` apskaičiavimas išskiriamas į `critical` bloką, tam, kad tik viena gija galėtų modifikuoti `X` reikšmę vienu metu. Apskaičiavus ir atnaujinus `X`, kiekviena gija susikuria savo `X` kopiją - `localX`. Šią kopiją galima naudoti `evaluateSolution` nes jinai ne bus keičiama. Kiekviena gija taip pat gauna `u` kopiją į kurią įrašo `evaluateSolution` apskaičiuotą reikšmę. Ciklo gale vėl naudojamas `critical`, tam kad tik viena gija vienu metu galėtų įvertinti `u > bestU` ir pakeisti `bestU` ir `bestX` reikšmes.

#figure(
  placement: none,
```c
bool increased = true;
int *manyXs = new int[NUM_THREADS * numX];

#pragma omp parallel private(u)
{
    while (increased) {
        int thread_id = omp_get_thread_num();
        int *localX = manyXs + (thread_id * numX);
        
        #pragma omp critical(increaseX)
        {
            increased = increaseX(X, numX-1, numCL);
            memcpy(localX, X, sizeof(int) * numX);
        }
        
        u = evaluateSolution(localX);

        #pragma omp critical(best)
        {
            if (u > bestU) {
                bestU = u;
                memcpy(bestX, localX, sizeof(int) * numX);
            }
        }
    }
}
```,
  caption: [Paralelizuotas visų galimų sprendinių perrinkimas]
)


=== Rezultatai

// === NUM_THREADS: 1 ===
#let values1 = (17.2353, 17.0992, 17.4169, 17.5637, 17.3402)
#let core1 = (values1.slice(0, count: 3).sum() / 3)

#let all_values1 = (21.8994, 21.8072, 22.1216, 22.2638, 21.9891)
#let all1 = (all_values1.slice(0, count: 3).sum() / 3)


// ```
// Matrix: 4.66416
// Spr: 17.2353
// All: 21.8994

// Matrix: 4.70804
// Spr: 17.0992
// All: 21.8072

// Matrix: 4.70465
// Spr: 17.4169
// All: 22.1216

// Matrix: 4.70014
// Spr: 17.5637
// All: 22.2638

// Matrix: 4.6489
// Spr: 17.3402
// All: 21.9891
// ```

// === NUM_THREADS: 2 ===
#let values2 = (9.29227, 9.08357, 9.14368, 9.09421, 9.13264)
#let core2 = (values2.slice(0, count: 3).sum() / 3)

#let all_values2 = (13.8621, 13.6069, 13.756, 13.4887, 13.7164)
#let all2 = (all_values2.slice(0, count: 3).sum() / 3)

// Matricos skaiciavimo trukme: 4.56983
// Sprendinio paieskos trukme: 9.29227
// Algoritmo vykdymo trukme: 13.8621

// Matricos skaiciavimo trukme: 4.52329
// Sprendinio paieskos trukme: 9.08357
// Algoritmo vykdymo trukme: 13.6069

// Matricos skaiciavimo trukme: 4.61234
// Sprendinio paieskos trukme: 9.14368
// Algoritmo vykdymo trukme: 13.756

// Matricos skaiciavimo trukme: 4.39453
// Sprendinio paieskos trukme: 9.09421
// Algoritmo vykdymo trukme: 13.4887

// Matricos skaiciavimo trukme: 4.58379
// Sprendinio paieskos trukme: 9.13264
// Algoritmo vykdymo trukme: 13.7164


// === NUM_THREADS: 4 ===
#let values4 = (4.77772, 4.73582, 4.85159, 4.63965, 4.51038)
#let core4 = (values4.slice(0, count: 3).sum() / 3)

#let all_values4 = (9.35896, 9.23228, 9.38276, 8.85537, 8.23031)
#let all4 = (all_values4.slice(0, count: 3).sum() / 3)

// Matricos skaiciavimo trukme: 4.58124
// Sprendinio paieskos trukme: 4.77772
// Algoritmo vykdymo trukme: 9.35896

// Matricos skaiciavimo trukme: 4.49646
// Sprendinio paieskos trukme: 4.73582
// Algoritmo vykdymo trukme: 9.23228

// Matricos skaiciavimo trukme: 4.53116
// Sprendinio paieskos trukme: 4.85159
// Algoritmo vykdymo trukme: 9.38276

// Matricos skaiciavimo trukme: 4.21572
// Sprendinio paieskos trukme: 4.63965
// Algoritmo vykdymo trukme: 8.85537

// Matricos skaiciavimo trukme: 3.71994
// Sprendinio paieskos trukme: 4.51038
// Algoritmo vykdymo trukme: 8.23031

#let beta = (core1 / all1)
#let alpha = 1 - beta

// Create the data for the two plots to overlay
#let data_core = (
  (1, core1 / core1), (2, core1 / core2), (4, core1 / core4)
)
#let data_all = (
  (1, all1 / all1), (2, all1 / all2), (4, all1 / all4)
)
#let data_linear = (
  (1, 1), (2, 2), (4, 4)
)
#let data_S_p = ((1, 1/(alpha + beta/1)), (2, 1/(alpha + beta/2)), (4, 1/(alpha + beta/4)))

// Create the axes for the overlay plot
#let x_axis = axis(min: 0, max: 4, step: 1, location: "bottom", title: "procesorių skaičius")
#let y_axis = axis(min: 0, max: 5, step: 1, location: "left", helper_lines: false, title: "pagreitėjimas")

// create a plot for each individual plot type and save the render call
#let pl_core = plot(data: data_core, axes: (x_axis, y_axis))
#let dp_core = graph_plot(pl_core, (50%, 25%), caption: "pagreitėjimo ir procesorių skaičiaus santykis", stroke: (paint: red, thickness: 1pt, dash: "dashed"))

#let pl_all = plot(data: data_all, axes: (x_axis, y_axis))
#let dp_all = graph_plot(pl_all, (50%, 25%), markings: "circle", stroke: (paint: green, thickness: 1pt, dash: "dashed"))

#let pl_linear = plot(data: data_linear, axes: (x_axis, y_axis))
#let dp_linear = graph_plot(pl_linear, (50%, 25%), stroke: blue)

#let pl_S_p = plot(data: data_S_p, axes: (x_axis, y_axis))
#let dp_S_p = graph_plot(pl_S_p, (50%, 25%), stroke: black)

#align(center)[#box()[
  #overlay((dp_core, dp_linear, dp_all, dp_S_p), (50%, 25%))
  #place(left, dx: 105%, dy: -87%)[Tiesinis pagreitėjimas]
  #place(left, dx: 105%, dy: -79%)[Sprendimo paieškos pagreitėjimas]
  // #place(left, dx: 105%, dy: -39.5%)[#mitex(`$\tilde{S}_p$`)]
  #place(left, dx: 105%, dy: -64%)[Teorinis pagreitėjimas]
  #place(left, dx: 105%, dy: -58%)[Programos pagreitėjimas]
]]

#pagebreak()
== Antra dalis

Šioje vietoje `for` direktyva atrodo lengvai pritaikoma, kadangi atstumų matricos kiekvieną eilutę galima apskaičiuoti nepriklausomai nuo to ar praeitos eilutės yra apskaičiuotos. Tačiau, pirmos eilutėms reikia žymiau mažiau laiko negu paskutinesnėms, todėl pritaikyta `guided` paskirstymo (angl. _scheduling_) direktyva. Tai leidžia efektyviau paskirstyti ciklų darbą per skirtingas gijas.

#figure(
  placement: none,
```c
distanceMatrix = new double*[numDP];
#pragma omp parallel
{
    #pragma omp for schedule(guided)
    for (int i=0; i<numDP; i++) {
        distanceMatrix[i] = new double[i+1];
        for (int j=0; j<=i; j++) {
            distanceMatrix[i][j] = HaversineDistance(
              demandPoints[i][0],
              demandPoints[i][1],
              demandPoints[j][0],
              demandPoints[j][1]
            );
        }
    }
}
}
```,
  caption: [Paralelizuotas atstumų matricos skaičiavimas]
) <matrix>

=== Rezultatai

/// MATRIX PARALLELIZATION

// === NUM_THREADS: 1 ===

#let values1 = (17.3815, 17.4154, 17.3834, 16.0638, 17.4336)
#let core1 = (values1.slice(0, count: 3).sum() / 3)

#let matrix_values1 = (4.53539, 4.54804, 4.54235, 4.59599, 4.54931)
#let matrix1 = (matrix_values1.slice(0, count: 3).sum() / 3)

#let all_values1 = (21.9169, 21.9635, 21.9257, 20.6598, 21.9829)
#let all1 = (all_values1.slice(0, count: 3).sum() / 3)

// Matricos skaiciavimo trukme: 4.53539
// Sprendinio paieskos trukme: 17.3815
// Algoritmo vykdymo trukme: 21.9169
// 
// Matricos skaiciavimo trukme: 4.54804
// Sprendinio paieskos trukme: 17.4154
// Algoritmo vykdymo trukme: 21.9635
// 
// Matricos skaiciavimo trukme: 4.54235
// Sprendinio paieskos trukme: 17.3834
// Algoritmo vykdymo trukme: 21.9257
// 
// Matricos skaiciavimo trukme: 4.59599
// Sprendinio paieskos trukme: 16.0638
// Algoritmo vykdymo trukme: 20.6598
// 
// Matricos skaiciavimo trukme: 4.54931
// Sprendinio paieskos trukme: 17.4336
// Algoritmo vykdymo trukme: 21.9829

// === NUM_THREADS: 2 ===

#let values2 = (8.71614, 9.21011, 9.12135, 8.66471, 8.94833)
#let core2 = (values2.slice(0, count: 3).sum() / 3)

#let matrix_values2 = (2.17127, 2.29965, 2.3158, 2.29705, 2.00836)
#let matrix2 = (matrix_values2.slice(0, count: 3).sum() / 3)

#let all_values2 = (10.8874, 11.5098, 11.4371, 10.9618, 10.9567)
#let all2 = (all_values2.slice(0, count: 3).sum() / 3)

// Matricos skaiciavimo trukme: 2.17127
// Sprendinio paieskos trukme: 8.71614
// Algoritmo vykdymo trukme: 10.8874
// 
// Matricos skaiciavimo trukme: 2.29965
// Sprendinio paieskos trukme: 9.21011
// Algoritmo vykdymo trukme: 11.5098
// 
// Matricos skaiciavimo trukme: 2.3158
// Sprendinio paieskos trukme: 9.12135
// Algoritmo vykdymo trukme: 11.4371
// 
// Matricos skaiciavimo trukme: 2.29705
// Sprendinio paieskos trukme: 8.66471
// Algoritmo vykdymo trukme: 10.9618
// 
// Matricos skaiciavimo trukme: 2.00836
// Sprendinio paieskos trukme: 8.94833
// Algoritmo vykdymo trukme: 10.9567

// === NUM_THREADS: 4 ===

#let values4 = (4.71473, 4.87095, 4.70493, 4.65552, 4.50287)
#let core4 = (values4.slice(0, count: 3).sum() / 3)

#let matrix_values4 = (1.08591, 1.04877, 1.03318, 1.05278, 1.06363)
#let matrix4 = (matrix_values4.slice(0, count: 3).sum() / 3)

#let all_values4 = (5.80063, 5.91972, 5.73811, 5.70831, 5.5665)
#let all4 = (all_values4.slice(0, count: 3).sum() / 3)

// Matricos skaiciavimo trukme: 1.08591
// Sprendinio paieskos trukme: 4.71473
// Algoritmo vykdymo trukme: 5.80063
// 
// Matricos skaiciavimo trukme: 1.04877
// Sprendinio paieskos trukme: 4.87095
// Algoritmo vykdymo trukme: 5.91972
// 
// Matricos skaiciavimo trukme: 1.03318
// Sprendinio paieskos trukme: 4.70493
// Algoritmo vykdymo trukme: 5.73811
// 
// Matricos skaiciavimo trukme: 1.05278
// Sprendinio paieskos trukme: 4.65552
// Algoritmo vykdymo trukme: 5.70831
// 
// Matricos skaiciavimo trukme: 1.06363
// Sprendinio paieskos trukme: 4.50287
// Algoritmo vykdymo trukme: 5.5665

#let beta = (core1 / all1)
#let alpha = 1 - beta

// Create the data for the two plots to overlay
#let data_core = (
  (1, core1 / core1), (2, core1 / core2), (4, core1 / core4)
)
#let data_all = (
  (1, all1 / all1), (2, all1 / all2), (4, all1 / all4)
)
#let data_linear = (
  (1, 1), (2, 2), (4, 4)
)
#let data_matrix = (
  (1, matrix1 / matrix1), (2, matrix1 / matrix2), (4, matrix1 / matrix4)
)

#let data_S_p = ((1, 1/(alpha + beta/1)), (2, 1/(alpha + beta/2)), (4, 1/(alpha + beta/4)))

// Create the axes for the overlay plot
#let x_axis = axis(min: 0, max: 4, step: 1, location: "bottom", title: "procesorių skaičius")
#let y_axis = axis(min: 0, max: 5, step: 1, location: "left", helper_lines: false, title: "pagreitėjimas")

// create a plot for each individual plot type and save the render call
#let pl_core = plot(data: data_core, axes: (x_axis, y_axis))
#let dp_core = graph_plot(pl_core, (50%, 25%), caption: "pagreitėjimo ir procesorių skaičiaus santykis", stroke: (paint: red, thickness: 1pt, dash: "dashed"))

#let pl_all = plot(data: data_all, axes: (x_axis, y_axis))
#let dp_all = graph_plot(pl_all, (50%, 25%), markings: "circle", stroke: (paint: green, thickness: 1pt, dash: "dashed"))

#let pl_linear = plot(data: data_linear, axes: (x_axis, y_axis))
#let dp_linear = graph_plot(pl_linear, (50%, 25%), stroke: blue)

#let pl_matrix = plot(data: data_matrix, axes: (x_axis, y_axis))
#let dp_matrix = graph_plot(pl_matrix, (50%, 25%), stroke: (paint: purple, thickness: 1pt, dash: "dashed"))

// #let pl_S_p = plot(data: data_S_p, axes: (x_axis, y_axis))
// #let dp_S_p = graph_plot(pl_S_p, (50%, 25%), stroke: black)

#align(center)[#box()[
  #overlay((dp_core, dp_linear, dp_all, dp_matrix), (50%, 25%))
  #place(left, dx: 100%, dy: -88%)[Tiesinis pagreitėjimas] // dp_linear
  #place(left, dx: 100%, dy: -76%)[Sprendimo paieškos pagreitėjimas] // do_core
  #place(left, dx: 100%, dy: -83%)[Programos pagreitėjimas] // dp_all
  #place(left, dx: 100%, dy: -94%)[Matricos skaičiavimo pagreitėjimas] // dp_matrix
]]
