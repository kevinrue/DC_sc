my_ggpairs <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
    geom_smooth(method=lm, fill="blue", color="blue", alpha = 0.5, ...) +
    geom_point()
  p
}
