limitRange <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping, ...) +
    geom_point(...) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0.25, 0.75)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0.25, 0.75))
}

limitRangeDens <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping, ...) +
    geom_density(...) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0.25, 0.75)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

my_fn <- function(data, mapping, method="p", use="pairwise", ...){

  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)

  # calculate correlation
  corr <- cor(x, y, method=method, use=use)

  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue,
  # zero to white, and one to red
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                              "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE",
                              "#4393C3", "#2166AC", "#053061"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]

  ggally_cor(data = data, mapping = mapping, digits = 2,
             stars = FALSE, size = 3, colour = I("black"),...) +
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

df_all <- as.data.frame(l_m_resamp_alphas[[1]])
colnames(df_all) <- c("alpha[HH]", "alpha[HW]", "alpha[HB]", "alpha[HO]",
                      "alpha[WH]", "alpha[WW]", "alpha[WB]", "alpha[WO]",
                      "alpha[BH]", "alpha[BW]", "alpha[BB]", "alpha[BO]",
                      "alpha[OH]", "alpha[OW]", "alpha[OB]", "alpha[OO]")

p1 <- ggpairs(df_all, cardinality_threshold = 16,
              upper = list(continuous = my_fn),
              lower = list(continuous = limitRange),
              diag = list(continuous = limitRangeDens), labeller = label_parsed) +
  labs(caption = "Note: H - Hispanic, W - NH White, B - NH Black, O - NH Other") +
  theme(plot.caption = element_text(hjust=0),
        strip.text = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(angle = 0))
print(p1)
ggsave("results/smdm_plot.png", p1, width = 12, height = 10)
ggsave("results/smdm_plot.jpg", p1, width = 12, height = 10)
ggsave("results/smdm_plot.pdf", p1, width = 12, height = 10)
ggsave("results/smdm_plot.tiff", p1, width = 12, height = 10, compression = "lzw")

jpeg("results/smdm_plot.jpeg", height = 800, width = 1000)
print(p1)
dev.off()

png("results/smdm_plot.png", height = 800, width = 1000)
print(p1)
dev.off()
