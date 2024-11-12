#import "@preview/mitex:0.2.4": *
#import "@preview/cetz:0.3.1": canvas, draw
#import "@preview/cetz-plot:0.1.0": plot, chart
#import "@preview/codly:1.0.0"
#import "data.typ": *

== `Attempt #5`

#let core1 = read_data(file: "../lab2/results/original.tsv", column: 2)
#let all1 = read_data(file: "../lab2/results/original.tsv", column: 3)
#let matrix1 = read_data(file: "../lab2/results/original.tsv", column: 1)

#let beta = (core1 / all1)
#let alpha = 1 - beta

#let data_linear = (
  (1, 1), (2, 2), (4, 4)
)
#let data_S_p = ((1, 1/(alpha + beta/1)), (2, 1/(alpha + beta/2)), (4, 1/(alpha + beta/4)), (8, 1/(alpha + beta/8)))

// =============================================

#let core2 = read_data(file: "../lab2/results/matrix_2.tsv", column: 2)
#let all2 = read_data(file: "../lab2/results/matrix_2.tsv", column: 3)
#let matrix2 = read_data(file: "../lab2/results/matrix_2.tsv", column: 1)

#let core4 = read_data(file: "../lab2/results/matrix_4.tsv", column: 2)
#let all4 = read_data(file: "../lab2/results/matrix_4.tsv", column: 3)
#let matrix4 = read_data(file: "../lab2/results/matrix_4.tsv", column: 1)

#let core8 = read_data(file: "../lab2/results/matrix_8.tsv", column: 2)
#let all8 = read_data(file: "../lab2/results/matrix_8.tsv", column: 3)
#let matrix8 = read_data(file: "../lab2/results/matrix_8.tsv", column: 1)

#let data_core = (
  (1, 1), (2, core1 / core2), (4, core1 / core4), (8, core1 / core8)
)

#let data_all = (
  (1, 1), (2, all1 / all2), (4, all1 / all4), (8, all1 / all8)
)
#let data_linear = (
  (1, 1), (2, 2), (4, 4), (8, 8)
)
#let data_matrix = (
  (1, 1), (2, matrix1 / matrix2), (4, matrix1 / matrix4), (8, matrix1 / matrix8)
)


#align(center)[
#figure(
  placement: none,
[
  #canvas({
    import draw: *

    // Set-up a thin axis style
    set-style(axes: (stroke: .5pt, tick: (stroke: .5pt)),
              legend: (stroke: none, orientation: ttb, item: (spacing: .3), scale: 70%))

    plot.plot(
      x-min: 0.9, x-max: 8.1,
      y-min: 0.0, y-max: 15.1,
      size: (10, 5.5),
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
          (x) => x - 1,
          domain: (1, 8),
          style: (stroke: (paint: rgb(0, 0, 0))), 
          label: "Teisinis pagreitėjimas (-1)"
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
  caption: [Pagreitėjimo ir Procesorių skaičiaus santykis (`first_attempt`)]
) <fig_first_attemp_1>
] 


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
      x-min: 0.9, x-max: 8.1,
      y-min: 0.0, y-max: 15.1,
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
          (x) => x - 1,
          domain: (1, 8),
          style: (stroke: (paint: rgb(0, 0, 0))), 
          label: "Teisinis pagreitėjimas (-1)"
        )

        plot.add(
          data_S_p,
          style: (stroke: (paint: rgb("#ff8104"))), 
          label: "Teorinis pagreitėjimas"
        )
        
        let _data_S_p = data_S_p.map((tuple) => (tuple.at(0) + 1, tuple.at(1)))
        plot.add(
          _data_S_p,
          style: (stroke: (paint: rgb("#ff0004"))), 
          label: "Teorinis pagreitėjimas (-1)"
        )
      })
  })

],
  supplement: "Fig. ",
  caption: [Pagreitėjimo ir Procesorių skaičiaus santykis (`first_attempt`)]
) <fig_first_attemp_2>


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
      x-min: 0.9, x-max: 8.1,
      y-min: 0.9, y-max: 15.2,
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

})],
  supplement: "Fig. ",
  caption: [Pagreitėjimo ir Procesorių skaičiaus santykis]
) <fig_first_attemp_2>