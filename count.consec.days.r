library(dplyr)

runchecker <- function(data, days){
  data %>% arrange(date) %>%
    group_by(id) %>%
    mutate(diff = c(0, diff(date)),
           periodID = 1 + cumsum(diff > days)) %>%
    group_by(id) %>%
    summarise(days = last(date) - first(date))
}

