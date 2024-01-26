library(tidyverse)
library(ggbeeswarm)

p1 <- read_csv("00_eQTL/results/model_AIC_dominant_recessive.csv") %>% 
  select(gene,snp, AIC_dominant, AIC_recessive) %>% 
  pivot_longer(-c(gene,snp)) %>% 
  mutate(name = recode(name,
                       "AIC_dominant"="Dominant",
                       "AIC_recessive"="Recessive")) %>% 
  
  ggplot(aes(x=name, y=value)) +
  geom_beeswarm(cex=3) +
  # geom_quasirandom(method = "smiley") +
  theme_classic(base_size = 10) +
  geom_hline(yintercept = c(-7,7), lty="dashed") +
  lims(y=c(-60,60)) +
  labs(x="Alternate model", 
       y="Change in AIC\nBetter fit by additive <---   ---> Better fit by alternate")
# p1

ggsave("publication/FigX_AIC.png", p1, width=2.5, height=4)
ggsave("publication/FigX_AIC.pdf", p1, width=2.5, height=4)
