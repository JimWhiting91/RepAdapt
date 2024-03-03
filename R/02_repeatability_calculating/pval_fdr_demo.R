# Script demonstrating different assumptions for p-values and fdr q-values
lib = c('data.table','qvalue','ggplot2','dplyr','tidyr')
sapply(lib,library,character.only = T)

# Draw either 1000,10000, or 100000 p-values from a random uniform distribution...
pval_size = c(1e3,1e4,1e5)
# Draw 1000 iterations of each
set.seed(12345)

all_sim_data = lapply(1:1000,function(i){
  pvals = lapply(pval_size,function(x){
    # Random uniform
    ps = runif(x)
    # FDR correct
    qs = p.adjust(ps,method = 'fdr')
    
    # How many p-values below <0.05
    signif_p = length(ps[ps < 0.05])
    
    # How many q-values below <0.5
    signif_fdr = length(qs[qs < 0.5])
    
    out = data.table(size = x,
                     `Signif p` = signif_p,
                     `Signif FDR` = signif_fdr,
                     iter = i)
  }) |>
    rbindlist()
}) |>
  rbindlist()

# Plot these out melted
all_sim_data$size = as.character(as.integer(all_sim_data$size))
all_sim_data = all_sim_data |>
  select(-iter) |>
  melt()

all_avg = all_sim_data[,.(median_signif = median(value)),by = .(size,variable)]

# Plot
demo_plot = ggplot(all_sim_data,aes(x = value)) +
  geom_histogram() +
  facet_grid(size~variable,scales = 'free') +
  geom_vline(data = all_avg,aes(xintercept = median_signif),colour = 'red2') +
  labs(x = 'Expected significant hits at p < 0.05 or fdr < 0.5',y = 'Count') +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 14))

# Save
pdf('figs/FigureSX_pval_fdr_demo.pdf',width = 6,height = 8)
demo_plot
dev.off()
