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

== Pagrindinės įžvalgos

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

Šį ciklą iš esmės yra ne paprasta parelilizuoti, nes:

  - `increaseX` - keičia masyvo `X` reikšmę (@increseX), ir ne tik `index`-ąjį elementą, rekursyviai kviesdamas save sumažina `index` reikšmę vienu, t.y. iškvietus `increaseX` visos masyvos reikšmės yra keičiamos. To pasekmė, kad `index`-ojo elemento skaičiavimo negalima paskirstyti skirtingoms gijoms, kitaip vėlesnėms gijoms reikėtų laukti, kol praeita gija baigs savo darbą, visiškai nustelbiant parelelizavimo naudą.

  - `increaseX` skaičiavimai priklauso vienas nuo kito, t.y. norint apskaičiuoti `X` reikšmę $n$-ame ciklo vykdyme, reikia pirma apskaičiuoti `X` reikšmę $(n-1)$-ame ciklo vykdyme. Analogiškai negalima paralelizuoti ir nes kitos gijos lauktų, kol praeita gija baigs savo darbą.

  - `if (u > bestU) { ... }` irgi gali tik vienas ciklas vienu metu, nes `bestU` ir `X` pakeitimas turi būti atliekamas "žingsniu".

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
  caption: [funkcijos `increseX` apibrėžimas]
) <increseX>

Tuo tarpu funkcija `evaluateSolution` nekeičia jokių globalių kintamųjų ar savo argumentų. Analitiškai žiūrint galima spėti, kad čia ir didžioji dalis skaičiavimo laiko yra sugaištama. Teorinis šios funkcijos "big-O" yra $O(#raw("numDP") dot #raw("numX"))$, tuo tarpu `increseX` rekursyviai save gali iškviesti daugiausiai `numX` kartų.

#figure(
  placement: none,
```c
double evaluateSolution(int *X) {
    double U = 0;
    double totalU = 0;
    int bestPF;
    int bestX;
    double d;

    #pragma omp parallel reduction (+:U, totalU) private(bestPF, bestX, d)
    #pragma omp for
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
  caption: [Paralelizuota funkcija `evaluateSolution`]
) <evaluateSolution>

@evaluateSolution esantį `for` ciklą galima parelilizuoti, kadangi `bestPF`, `bestX`, `d` kintamiesiems ciklo viduje prieš jų panaudojimą priskiriamos tų pačių konstantų reikšmės ir kaip minėta `X` reikšmė nekeičiama. Kiekvienai gijai sukurianant atskirą `bestPF`, `bestX`, `d` kopiją išvengiamos "data-race" problemos (panaudojant `private` direktyvą). Už `for` ciklo ribų reikšmingi tik `totalU` ir `U` kintamieji, juos apsaugi nuo "data-race" galima apsaugti panaudojant `reduce` direktyvą, kadangi jie naudojami tik galutiniui rezultatui susumuoti.

== Rezultatai \#1

// === NUM_THREADS: 1 ===
#let values1 = (17.2353, 17.0992, 17.4169, 17.5637, 17.3402)
#let core1 = (values1.slice(0, count: 3).sum() / 3)

// NUM_THREADS = 1

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


// Create the data for the two plots to overlay
#let data_scatter = (
  (1, core1 / core1), (2, core1 / core2), (4, core1 / core4)
)
#let data_graph = (
  (0, 1), (1, 2), (2, 1)
)

// Create the axes for the overlay plot
#let x_axis = axis(min: 0, max: 4, step: 1, location: "bottom")
#let y_axis = axis(min: 0, max: 5, step: 1, location: "left", helper_lines: false)

// create a plot for each individual plot type and save the render call
#let pl_scatter = plot(data: data_scatter, axes: (x_axis, y_axis))
#let scatter_display = scatter_plot(pl_scatter, (50%, 25%), stroke: red)
#let pl_graph = plot(data: data_graph, axes: (x_axis, y_axis))
#let graph_display = graph_plot(pl_graph, (50%, 25%), stroke: blue)

// overlay the plots using the overlay function
#overlay((scatter_display, graph_display), (50%, 25%))

/// MATRIX PARALLELIZATION

// === NUM_THREADS: 1 ===
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
// 
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

// === NUM_THREADS: 2 ===
//
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
