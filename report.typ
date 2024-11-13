#import "@preview/mitex:0.2.4": *
#import "@preview/cetz:0.3.1": canvas, draw
#import "@preview/cetz-plot:0.1.0": plot, chart
#import "@preview/codly:1.0.0"
#import "typst/data.typ": *

#show: codly.codly-init.with()
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

#align(center, text(26pt)[
  *Laboratoriniai darbai*
])

#align(center, text(15pt)[
  Domantas Keturakis \
  #text(13pt)[
    Lapkritis 2024
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

#show outline.entry.where(
  level: 1
): it => {
  v(1em, weak: true)
  strong(underline(smallcaps(it.body)))
}

#show outline.entry.where(
  level: 2
): it => {
  strong(it.body)
}

#show outline.entry.where(
  level: 4
): it => {
  v(0.5em, weak: true)
  emph(it.body)
}

#outline(
  indent: 0.5em,
  title: [Turinys],
)

#pagebreak()

= Užduotis \#1

== Pirma dalis

=== Lygiagretinimo galimybių analizė

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

  - `increaseX` skaičiavimai priklauso vienas nuo kito, t.y. norint apskaičiuoti `X` reikšmę $n$-ame ciklo vykdyme, reikia pirma apskaičiuoti `X` reikšmę $(n-1)$-ame ciklo vykdyme. Analogiškai negalima lygiagretinti nes kitos gijos lauktų, kol praeita gija baigs savo darbą.

  - `if (u > bestU) { ... }` irgi gali tik vienas ciklas vienu metu, nes `bestU` ir `X` pakeitimas turi būti atliekamas "žingsniu" - t.y. atomiškai.

Iš esmės neperrašius `increaseX`, `X` skaičiavimų lygiagretinti nėra praktiška.

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
    double U = 0; double totalU = 0;
    int bestPF, bestX;
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


=== Sprendimo lygiagretinimas

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
  caption: [Sulygiagretintas visų galimų sprendinių perrinkimas]
)

=== Rezultatai

// #let results = csv("lab1/no_matrix_1.tsv", delimiter: "\t")
// #table(
//   columns: 4,
//   rows: 4,
//   [*Threads*],
//   [$"t_matrix" - "t_start"$],
//   [$"t_finish" - "t_matrix"$],
//   [$"t_finish" - "t_start"$],
//   ..results.flatten(),
// )

#let core1 = read_data(file: "../lab1/no_matrix_1.tsv", column: 2)
#let all1 = read_data(file: "../lab1/no_matrix_1.tsv", column: 3)

#let core2 = read_data(file: "../lab1/no_matrix_2.tsv", column: 2)
#let all2 = read_data(file: "../lab1/no_matrix_2.tsv", column: 3)

#let core4 = read_data(file: "../lab1/no_matrix_4.tsv", column: 2)
#let all4 = read_data(file: "../lab1/no_matrix_4.tsv", column: 3)


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
#let data_S_p = ((1, 1/(alpha + beta/1)), (2, 1/(alpha + beta/2)), (4, 1/(alpha + beta/4)), (8, 1/(alpha + beta/8)))

#figure(
  placement: none,
[
  #canvas({
    import draw: *

    // Set-up a thin axis style
    set-style(
      axes: (stroke: .5pt, tick: (stroke: .5pt)),
      legend: (stroke: none, orientation: ttb, item: (spacing: .3), scale: 80%),
    )

    plot.plot(
      x-min: 0.9, x-max: 4.2,
      y-min: 0.9, y-max: 4.5,
      size: (10, 6),
      x-tick-step: 1,
      y-tick-step: 1,
      y-minor-tick-step: 0.5,
      x-label: [Procesorių skaičius],
      y-label: [Pagreitėjimas],
      x-grid: true,
      y-grid: true,
      // legend: auto,
      {

        plot.add(
          data_core,
          style: (stroke: (paint: green, dash: "dashed")), 
          label: "Sprendimo paieškos pagreitėjimas"
        )

        plot.add(
          data_all,
          style: (stroke: (paint: rgb("#e64914"), dash: "dashed")), 
          label: "Programos pagreitėjimas"
        )

        plot.add(
          data_linear,
          style: (stroke: (paint: blue)), 
          label: "Tiesinis pagreitėjimas"
        )

        plot.add(
          data_S_p,
          style: (stroke: (paint: rgb("#ff8104"))), 
          label: "Teorinis pagreitėjimas"
        )
      })
  })

],
  supplement: "Fig. ",
  caption: [Pagreitėjimo ir Procesorių skaičiaus santykis, kai matricos skaičiavimas nuoseklus]
) <fig1>

Iš @fig1 matoma, kad lygiagretinamos dalies pagreitėjimas seka tiesinį pagreitėjimą. O visos programos pagreitėmijas seka teorinį pagreitėjimą nusakytą pagal Amdalo dėsnį.

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
```,
  caption: [Sulygiagretintas atstumų matricos skaičiavimas]
) <matrix>

=== Rezultatai

#let core1 = read_data(file: "../lab1/with_matrix_1.tsv", column: 2)
#let all1 = read_data(file: "../lab1/with_matrix_1.tsv", column: 3)
#let matrix1 = read_data(file: "../lab1/with_matrix_1.tsv", column: 1)

#let core2 = read_data(file: "../lab1/with_matrix_2.tsv", column: 2)
#let all2 = read_data(file: "../lab1/with_matrix_2.tsv", column: 3)
#let matrix2 = read_data(file: "../lab1/with_matrix_2.tsv", column: 1)

#let core4 = read_data(file: "../lab1/with_matrix_4.tsv", column: 2)
#let all4 = read_data(file: "../lab1/with_matrix_4.tsv", column: 3)
#let matrix4 = read_data(file: "../lab1/with_matrix_4.tsv", column: 1)

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


#align(center)[
#figure(
  placement: none,
[
  #canvas({
    import draw: *

    // Set-up a thin axis style
    set-style(axes: (stroke: .5pt, tick: (stroke: .5pt)),
              legend: (stroke: none, orientation: ttb, item: (spacing: .3), scale: 80%))

    plot.plot(
      x-min: 0.9, x-max: 4.1,
      y-min: 0.9, y-max: 4.1,
      size: (10, 6),
      x-tick-step: 1,
      y-tick-step: 1,
      y-minor-tick-step: 0.5,
      x-label: [Procesorių skaičius],
      y-label: [Pagreitėjimas],
      x-grid: true,
      y-grid: true,
      {

        plot.add(
          data_core,
          style: (stroke: (paint: green, dash: "dashed")), 
          label: "Sprendimo paieškos pagreitėjimas"
        )

        plot.add(
          data_all,
          style: (stroke: (paint: rgb("#e64914"), dash: "dashed")), 
          label: "Programos pagreitėjimas"
        )

        plot.add(
          data_linear,
          style: (stroke: (paint: blue)), 
          label: "Tiesinis pagreitėjimas"
        )

        plot.add(
          data_matrix,
          // mark: "x", mark-size: 0.15,
          style: (stroke: (paint: orange)), 
          label: "Matricos skaičiavimo pagreitėjimas",
        )
      })
  })
  
],
  supplement: "Fig. ",
  caption: [Pagreitėjimo ir Procesorių skaičiaus santykis, kai matricos skaičiavimas sulygiagretintas]
) <fig2>
] 

Iš @fig2 matoma, kad matricos skaičiavimo ir sprendimo ieškojimo pagreitėjimas seka tiesinę kreivę, kas matosi ir visos programos pagreitėjime, kuri irgi seka tiesinę kreivę.

#pagebreak()
= Užduotis \#2

Sulygiagretinti skaičiavimus naudojantis `MPI` biblioteka.

== Pirma dalis

Sulygiagretinti sprendinio skaičiavimus.

=== Pirmas bandymas

Procesus galima paskirstyti taip, kad yra vienas pagrindinis (_main_) procesas ir likę - darbuotojai (_workers_).
Kur pagrindinis procesas skaičiuoja atnaujintas `X` reikšmes ir jas išsiunčia kitiem procesam.

```c
int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

// Iškirpta: loadDemandPoints(), distanceMatrix skaičiavimas, MPI_Buffer_attach, t.t.
MPI_Barrier(MPI_COMM_WORLD);

bool increased = true;
while (increased) {
   int sends = 1;
   if (world_rank == 0) {
      for(int ix = 1; ix < world_size; ++ix) {
         increased = increaseX(X, numX - 1, numCL);

         MPI_Send(X, numX, MPI_INT, ix, DATA_X, MPI_COMM_WORLD);
         sends = ix + 1;

         if (!increased) {
            for (int ix = 1; ix < world_size; ++ix) {
               MPI_Bsend(X, 0, MPI_BYTE, ix, SIGNAL_DONE, MPI_COMM_WORLD);
            }
            break;
         }
      }
   }
   // ...
}
```

Radus galutinę `X` reikšmę pagrindinis procesas išsiučia žinutę su žyme `SIGNAL_DONE`, kad pranešti procesams darbuotojiems, kad šie gali baigti savo darbą.

#pagebreak()

Tuo tarpu procesai darbuotojai laukia naujos `X` reiškės, jos sulaukę, apskaičiuoja `u` reikšmę ir išsiunčia ją, kartu su `X`, pagrindiniam procesui. Taip pat jie patikrinina ar negavo žinutės su `SIGNAL_DONE` žyma, kurios pagalba jie sužino, kad daugiau negaus `X` reiškmių ir todėl gali baigti savo darbą.

```c
while(increased) {
  // ...
  if (world_rank != 0) {
    MPI_Status status;
    MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (status.MPI_TAG == SIGNAL_DONE) {
     MPI_Recv(NULL, 0, MPI_BYTE, 0, SIGNAL_DONE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     break;
    }

    MPI_Recv(X, numX, MPI_INT, 0, DATA_X, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    u = evaluateSolution(X);

    MPI_Send(&u, 1, MPI_DOUBLE, 0, DATA_U, MPI_COMM_WORLD);
    MPI_Send(X, numX, MPI_INT,  0, DATA_X, MPI_COMM_WORLD);
  }
  // ...
}
```

Tuo tarpu pagrindinis laukia tiek `u` ir `X` reikšmių kiek išsiuntė procesams darbuotojams, gavęs reikšmę su didesne `u` reikšme atnauja savo `bestU` ir `bestX` kintamuosius.

```c
while(increased) {
   // ...
   if (world_rank == 0) {
      for (int ix = 1; ix < sends; ++ix) {
         MPI_Recv(&u, 1, MPI_DOUBLE, ix, DATA_U, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         MPI_Recv(X, numX, MPI_INT, ix, DATA_X, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

         if (u > bestU) {
            bestU = u;
            memcpy(bestX, X, sizeof(int) * numX);
         }
      }
   }
}
```

==== Rezultatai

#let core1 = read_data(file: "../lab2/results/1_original.tsv", column: 2)
#let all1 = read_data(file: "../lab2/results/1_original.tsv", column: 3)
#let matrix1 = read_data(file: "../lab2/results/1_original.tsv", column: 1)

#let beta = (core1 / all1)
#let alpha = 1 - beta
#let data_S_p = ((1, 1/(alpha + beta/1)), (2, 1/(alpha + beta/2)), (4, 1/(alpha + beta/4)), (8, 1/(alpha + beta/8)))

#let data_linear = (
  (1, 1), (2, 2), (4, 4), (8, 8)
)

// =============================================

#let core2 = read_data(file: "../lab2/results/1_2.tsv", column: 2)
#let all2 = read_data(file: "../lab2/results/1_2.tsv", column: 3)
#let matrix2 = read_data(file: "../lab2/results/1_2.tsv", column: 1)

#let core4 = read_data(file: "../lab2/results/1_4.tsv", column: 2)
#let all4 = read_data(file: "../lab2/results/1_4.tsv", column: 3)
#let matrix4 = read_data(file: "../lab2/results/1_4.tsv", column: 1)

#let core8 = read_data(file: "../lab2/results/1_8.tsv", column: 2)
#let all8 = read_data(file: "../lab2/results/1_8.tsv", column: 3)
#let matrix8 = read_data(file: "../lab2/results/1_8.tsv", column: 1)

#let data_core = (
  (2, core1 / core2), (4, core1 / core4), (8, core1 / core8)
)
#let data_all = (
  (2, all1 / all2), (4, all1 / all4), (8, all1 / all8)
)
#let data_matrix = (
  (2, matrix1 / matrix2), (4, matrix1 / matrix4), (8, matrix1 / matrix8)
)

#grid(
  columns: 2,
  gutter: 2mm,
[#figure(
  placement: none,
[
  #canvas({
    import draw: *

    // Set-up a thin axis style
    set-style(
      axes: (stroke: .5pt, tick: (stroke: .5pt)),
      legend: (stroke: none, fill: none, orientation: ttb, item: (spacing: .1), scale: 40%),
    )

    plot.plot(
      x-min: 0.9, x-max: 8.1,
      y-min: 0.9, y-max: 8.1,
      size: (7, 5),
      x-tick-step: 1,
      y-tick-step: 1,
      y-minor-tick-step: 0.5,
      x-label: [Procesorių skaičius],
      y-label: [Pagreitėjimas],
      x-grid: true,
      y-grid: true,
      legend: "inner-north-west",
      {

        plot.add(
          data_core,
          style: (stroke: (paint: green, dash: "dashed")), 
          label: text(8pt)[Sprendimo paieškos pagreitėjimas],
        )

        plot.add(
          data_all,
          style: (stroke: (paint: rgb("#e64914"), dash: "dashed")), 
          label: text(8pt)[Programos pagreitėjimas],
        )

        plot.add(
          data_linear,
          style: (stroke: (paint: blue)), 
          label: text(8pt)[Tiesinis pagreitėjimas],
        )

        plot.add(
          data_S_p,
          style: (stroke: (paint: rgb("#ff8104"))), 
          label: text(8pt)[Teorinis pagreitėjimas],
        )
      })
  })

],
  supplement: "Fig. ",
  caption: [Pagreitėjimo ir Procesorių skaičiaus santykis]
)<lab2_fig_1>],
[
#let core3 = read_data(file: "../lab2/results/1_3.tsv", column: 2)
#let all3 = read_data(file: "../lab2/results/1_3.tsv", column: 3)
#let matrix3 = read_data(file: "../lab2/results/1_3.tsv", column: 1)

#let core5 = read_data(file: "../lab2/results/1_5.tsv", column: 2)
#let all5 = read_data(file: "../lab2/results/1_5.tsv", column: 3)
#let matrix5 = read_data(file: "../lab2/results/1_5.tsv", column: 1)

#let core9 = read_data(file: "../lab2/results/1_9.tsv", column: 2)
#let all9 = read_data(file: "../lab2/results/1_9.tsv", column: 3)
#let matrix9 = read_data(file: "../lab2/results/1_9.tsv", column: 1)

#let data_core = (
  (2, core1 / core3), (4, core1 / core5), (8, core1 / core9)
)

#let data_all = (
  (2, all1 / all3), (4, all1 / all5), (8, all1 / all9)
)
#let data_matrix = (
  (2, matrix1 / matrix3), (4, matrix1 / matrix5), (8, matrix1 / matrix9)
)

#figure(
  placement: none,
[
  #canvas({
    import draw: *

    // Set-up a thin axis style
    set-style(
      axes: (stroke: .5pt, tick: (stroke: .5pt)),
      legend: (stroke: none, fill: none, orientation: ttb, item: (spacing: .1), scale: 40%),
    )

    plot.plot(
      x-min: 0.9, x-max: 8.1,
      y-min: 0.9, y-max: 8.1,
      size: (7, 5),
      x-tick-step: 1,
      y-tick-step: 1,
      y-minor-tick-step: 0.5,
      x-label: [Procesorių skaičius],
      y-label: [Pagreitėjimas],
      x-grid: true,
      y-grid: true,
      legend: "inner-north-west",
      {

        plot.add(
          data_core,
          style: (stroke: (paint: green, dash: "dashed")), 
          label: text(8pt)[Sprendimo paieškos pagreitėjimas]
        )

        plot.add(
          data_all,
          style: (stroke: (paint: rgb("#e64914"), dash: "dashed")), 
          label: text(8pt)[Programos pagreitėjimas]
        )

        plot.add(
          (x) => x,
          domain: (1, 8),
          style: (stroke: (paint: blue)), 
          label: text(8pt)[Teisinis pagreitėjimas (-1)]
        )
        
        plot.add(
          data_S_p,
          style: (stroke: (paint: rgb("#ff8104"))), 
          label: text(8pt)[Teorinis pagreitėjimas (-1)]
        )
      })
  })

],
  supplement: "Fig. ",
  caption: [Pagreitėjimo ir Procesorių skaičiaus santykis, neskaičiuojant pagrindinio proceso]
) <lab2_fig_1_without_main>
])

Šis sprendimas nėra prastas (@lab2_fig_1), bet pilnai neišnaudoja pagrindinio proceso, galima daryti išvadas, kad jis dažnai neturi darbo ir laukio kol galės kitiem procesam išsiųsti atnaujintas `X` reikšmes. Neskaičiuojant pagrindinio proceso (@lab2_fig_1_without_main), t.y. skaičiuojant pagreitėjimą su 8 procesais ištikrųjų paleidžiami 9 procesai, praktiškai pasiekiamas teorinis pagreitėjimas.

==== Antras bandymas

Visgi galima padaryti taip, kad pagrindinis procesas irgi skaičiuotu `X` reikšmes. Kad procesai darbuotojai nelauktų, kol pagrindinis procesas baigs skaičiuoti savo dalį, nebelaukiama atsakymo iš procesų darbuotojų, `X` siuntimui pasitelkiamas `MPI_Bsend`.

```c
while (true) {
   if (world_rank == 0 && increased) {
      for(int ix = 1; ix < world_size; ++ix) {
         increased = increaseX(X, numX - 1, numCL);
         MPI_Bsend(X, numX, MPI_INT, ix, DATA_X, MPI_COMM_WORLD);
         sends += 1;

         if (!increased) {
            for (int ix = 1; ix < world_size; ++ix) {
               MPI_Bsend(&dummy_load, 1, MPI_INT, ix, SIGNAL_DONE, MPI_COMM_WORLD);
            }
            break;
         }
      }
   }
// ...
}
```

Toliau `while` cikle pridedamas paskaičiavimas pagrindiniam procesui ir `SIGNAL_DONE` išsiuntimas, jeigu šiame cikle baigtusi `X` skaičiavimas:

```c
while(true) {
   // ...
   if (world_rank == 0 && increased) {
      increased = increaseX(X, numX - 1, numCL);
      
      u = evaluateSolution(X);
      if (u > bestU) {
         bestU = u;
         memcpy(bestX, X, sizeof(int) * numX);
      }

      if (!increased) {
         for (int ix = 1; ix < world_size; ++ix) {
            MPI_Bsend(&dummy_load, 1, MPI_INT, ix, SIGNAL_DONE, MPI_COMM_WORLD);
         }
      }
   }
   // ...
}
```

Procesams darbuotojams pagrinde ne daug kas keičiasi. Jie irgi pakeičia `MPI_Send` į `MPI_Bsend`.

```c
   if (world_rank != 0) {
      int master_done = WRP_Check_for(0, SIGNAL_DONE, MPI_COMM_WORLD);
      int master_sent_X = WRP_Check_for(0, DATA_X, MPI_COMM_WORLD);

      if (master_sent_X) {
         MPI_Recv(X, numX, MPI_INT, 0, DATA_X, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

         u = evaluateSolution(X);
         if (bestU < u) { bestU = u; bestX = memcpy(bestX, X, sizeof(int) * numX); }

         MPI_Bsend(&u, 1, MPI_DOUBLE, 0, DATA_U, MPI_COMM_WORLD);
         MPI_Bsend(X, numX, MPI_INT,  0, DATA_X, MPI_COMM_WORLD);
      }

      if (master_done && !master_sent_X) { break; }
   }
```

Čia `WRP_Check_for(int source, int tag, MPI_Comm comm)` (apibrėžimas @WRP_Check_for) viduje naudoja `MPI_Iprobe(int source, int tag, MPI_Comm communicator, int* flag, MPI_Status* status)` ir gražina `flag` dalį.

```c
   if (world_rank == 0) {
      if (!increased && receives == sends) { break; }

      MPI_Status status;
      int worker_sent_X;
      MPI_Iprobe(MPI_ANY_SOURCE, DATA_X, MPI_COMM_WORLD, &worker_sent_X, &status);
      while(worker_sent_X) {
         MPI_Recv(&copy_u, 1, MPI_DOUBLE, status.MPI_STATUS, DATA_U, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         MPI_Recv(copy_X, numX, MPI_INT, status.MPI_STATUS, DATA_X, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         receives += 1;

         if (copy_u > bestU) {
            bestU = copy_u;
            memcpy(bestX, copy_X, sizeof(int) * numX);
         }

         worker_sent_X = WRP_Check_for(MPI_ANY_SOURCE, DATA_X, MPI_COMM_WORLD);
      }
   }
```


#let core2 = read_data(file: "../lab2/results/4_2.tsv", column: 2)
#let all2 = read_data(file: "../lab2/results/4_2.tsv", column: 3)
#let matrix2 = read_data(file: "../lab2/results/4_2.tsv", column: 1)

#let core4 = read_data(file: "../lab2/results/4_4.tsv", column: 2)
#let all4 = read_data(file: "../lab2/results/4_4.tsv", column: 3)
#let matrix4 = read_data(file: "../lab2/results/4_4.tsv", column: 1)

#let core8 = read_data(file: "../lab2/results/4_8.tsv", column: 2)
#let all8 = read_data(file: "../lab2/results/4_8.tsv", column: 3)
#let matrix8 = read_data(file: "../lab2/results/4_8.tsv", column: 1)

#let data_core = (
  (2, core1 / core2), (4, core1 / core4), (8, core1 / core8)
)
#let data_all = (
  (2, all1 / all2), (4, all1 / all4), (8, all1 / all8)
)
#let data_matrix = (
  (2, matrix1 / matrix2), (4, matrix1 / matrix4), (8, matrix1 / matrix8)
)

#figure(
  placement: none,
[
  #canvas({
    import draw: *

    // Set-up a thin axis style
    set-style(
      axes: (stroke: .5pt, tick: (stroke: .5pt)),
      legend: (stroke: none, fill: none, orientation: ttb, item: (spacing: .1), scale: 70%),
    )

    plot.plot(
      x-min: 0.9, x-max: 8.1,
      y-min: 0.9, y-max: 8.1,
      size: (10, 6),
      x-tick-step: 1,
      y-tick-step: 1,
      y-minor-tick-step: 0.5,
      x-label: [Procesorių skaičius],
      y-label: [Pagreitėjimas],
      x-grid: true,
      y-grid: true,
      legend: "inner-north-west",
      {

        plot.add(
          data_core,
          style: (stroke: (paint: green, dash: "dashed")), 
          label: text(8pt)[Sprendimo paieškos pagreitėjimas],
        )

        plot.add(
          data_all,
          style: (stroke: (paint: rgb("#e64914"), dash: "dashed")), 
          label: text(8pt)[Programos pagreitėjimas],
        )

        plot.add(
          data_linear,
          style: (stroke: (paint: blue)), 
          label: text(8pt)[Tiesinis pagreitėjimas],
        )

        plot.add(
          data_S_p,
          style: (stroke: (paint: rgb("#ff8104"))), 
          label: text(8pt)[Teorinis pagreitėjimas],
        )
      })
  })

],
  supplement: "Fig. ",
  caption: [Pagreitėjimo ir Procesorių skaičiaus santykis]
)<lab2_fig_4>

Šis metodas pasieka teorinį pagreitėjimą be papildomo proceso.

#pagebreak()
== Antra dalis

=== Užduotis

Sulygiagretinti matricos skaičiavimus.

=== Lygiagretinimas

Tam, kad gerai sulygiagretinti, reikia paskirstyti skaičiavimus taip, kad kiekvienas procesas tūrėtų po tiek pat darbo. Tam panaudojama funkcija `lengths` (Jos apibrėžimas - @choose_interval), kuri parenka tinkamus intervalus, taip, kad kiekvienas procesas paskaičiuotu apytiksliai tiek pat matricos elementų.

```c
int main() {
   // ...
   int *lens = lengths(numDP, world_size);
   int *counts = calloc(sizeof(int), world_size);
   for (int ix = 0; ix < world_size; ++ix) {
      counts[ix] = (lens[ix + 1] - lens[ix]) * numDP;
   }

   int *disps = calloc(sizeof(int), world_size);
   for (int ix = 0; ix < world_size; ++ix) {
      disps[ix] = lens[ix] * numDP;
   }

   distanceMatrix = calloc(sizeof(double), numDP * numDP);
   for (int i = lens[world_rank]; i < lens[world_rank + 1]; i++) {
      for (int j = 0; j <= i; j++) {
         distanceMatrix[numDP * i + j] =
            HaversineDistance4(demandPoints[i][0], demandPoints[i][1], demandPoints[j][0], demandPoints[j][1]);
      }
   }
   // ...
```

Kai procesas baigia savo dalį, jis nusiunčia nusiunčia kitiem procesam savo dalį ir laukia kitų dalių naudojant `MPI_Allgatherv`.

```c
   MPI_Allgatherv(
      distanceMatrix + disps[world_rank],
      counts[world_rank],
      MPI_DOUBLE,
      distanceMatrix,
      counts,
      disps,
      MPI_DOUBLE,
      MPI_COMM_WORLD
   )
```

#pagebreak()
=== Rezultatai

TODO: actual data

#let core2 = read_data(file: "../lab2/results/4_2.tsv", column: 2)
#let all2 = read_data(file: "../lab2/results/4_2.tsv", column: 3)
#let matrix2 = read_data(file: "../lab2/results/4_2.tsv", column: 1)

#let core4 = read_data(file: "../lab2/results/4_4.tsv", column: 2)
#let all4 = read_data(file: "../lab2/results/4_4.tsv", column: 3)
#let matrix4 = read_data(file: "../lab2/results/4_4.tsv", column: 1)

#let core8 = read_data(file: "../lab2/results/4_8.tsv", column: 2)
#let all8 = read_data(file: "../lab2/results/4_8.tsv", column: 3)
#let matrix8 = read_data(file: "../lab2/results/4_8.tsv", column: 1)

#let data_core = (
  (2, core1 / core2), (4, core1 / core4), (8, core1 / core8)
)
#let data_all = (
  (2, all1 / all2), (4, all1 / all4), (8, all1 / all8)
)
#let data_matrix = (
  (2, matrix1 / matrix2), (4, matrix1 / matrix4), (8, matrix1 / matrix8)
)

#align(center)[
#figure(
  placement: none,
[
  #canvas({
    import draw: *

    // Set-up a thin axis style
    set-style(axes: (stroke: .5pt, tick: (stroke: .5pt)),
              legend: (stroke: none, orientation: ttb, item: (spacing: .3), scale: 80%))

    plot.plot(
      x-min: 0.9, x-max: 8.1,
      y-min: 0.9, y-max: 8.1,
      size: (10, 6),
      x-tick-step: 1,
      y-tick-step: 1,
      y-minor-tick-step: 0.5,
      x-label: [Procesorių skaičius],
      y-label: [Pagreitėjimas],
      x-grid: true,
      y-grid: true,
      {

        plot.add(
          data_core,
          style: (stroke: (paint: green, dash: "dashed")), 
          label: "Sprendimo paieškos pagreitėjimas"
        )

        plot.add(
          data_all,
          style: (stroke: (paint: rgb("#e64914"), dash: "dashed")), 
          label: "Programos pagreitėjimas"
        )

        plot.add(
          (x) => x,
          domain: (1, 8),
          style: (stroke: (paint: blue)), 
          label: "Tiesinis pagreitėjimas"
        )

        plot.add(
          data_matrix,
          // mark: "x", mark-size: 0.15,
          style: (stroke: (paint: orange)), 
          label: "Matricos skaičiavimo pagreitėjimas",
        )
      })
  })
  
],
  supplement: "Fig. ",
  caption: [Pagreitėjimo ir Procesorių skaičiaus santykis, kai matricos skaičiavimas sulygiagretintas]
) <lab2_fig_matrix>
] 

#pagebreak()
= Priedai <appendixes>

#figure(
  placement: none,
```c
int WRP_Check_for(int source, int tag, MPI_Comm comm) {
   MPI_Status status;
   int flag;
   MPI_Iprobe(source, tag, comm, &flag, &status);

   return flag;
}
```,
  caption: [Optimalus intervalus parinkimas]
) <WRP_Check_for>

#figure(
  placement: none,
```c
int *lengths(int leg_length, int process_count) {
    int area = leg_length * leg_length / 2;

    int small_area = area / process_count;
    int *lenghts = calloc(sizeof(int), process_count + 1);
    for (int ix = 1; ix < process_count; ++ix) {
        lenghts[ix] = (int) sqrt(2 * ( ix ) * small_area);
    }
    lenghts[process_count] = leg_length;

    return lenghts;
}
```,
  caption: [Optimalus intervalus parinkimas]
) <choose_interval>


// #include "typst/2_attempt.typ"

// #include "typst/3_attempt.typ"

// #include "typst/4_attempt.typ"

// #include "typst/5_matrix.typ"

// #include "typst/6_matrix2.typ"

// #include "typst/6_matrix3.typ"
