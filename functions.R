set_panel_heights <- function(g, heights){
  g$heights <- grid:::unit.list(g$heights) # hack until R 3.3 comes out
  id_panels <- unique(g$layout[g$layout$name=="panel", "t"])
  g$heights[id_panels] <- heights
  g
}