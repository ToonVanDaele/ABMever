# Set hunting scenario's

get_hunting_scen <- function(file){

  path <- "./data/input/hunting_scenarios.xlsx"

  Hscen <- path %>%
    readxl::excel_sheets() %>%
    set_names() %>%
    map(readxl::read_excel, path = path) %>%
  return(Hscen)
}
