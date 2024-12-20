
#let read_data(file: "", column: 1, average_count: 3) = {
  let results = csv(file, delimiter: "\t")

  let values = results
    .map(r => float(r.at(column)))
    .slice(0, count: average_count)
    .sum()
  
  values / average_count
}

// #let data(files: (), column: 1) = {
//   let numbers = (1, 2, 3, 4, 5, 6,)
//   let records = files.map(f => read_data(file: f, column: column))

//   let base_record = records.at(0)
//   records.map(r => r / base_record)
// }