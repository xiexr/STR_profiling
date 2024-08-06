###################################

library(extrafont)
library(ggplot2)

# Clear the environment
rm(list = ls())

font_import(pattern = "arial", prompt = FALSE)
loadfonts(device = "win")
# Set the working directory
setwd('E:/experiment/Ğ£×¼´íÎóÂÊ/Nip')

# Read the data from the file
data_gatk  <- read.csv('gatkSTRrepeattimecount.txt', sep = '\t')
data_hip <- read.csv('hipSTRrepeattimecount.txt', sep = '\t')
data_lob <- read.csv('lobSTRrepeattimecount.txt', sep = '\t')

x_range <- range(c(data_gatk$refrptime, data_hip$refrptime, data_lob$refrptime))
y_range <- range(c(data_gatk$gatkSTRrepeattime, data_hip$hiprepeattime, data_lob$lobrepeattime))

# Create a heatmap with circles
scatter_plot <- ggplot(data_gatk, aes(x = refrptime, y = gatkSTRrepeattime)) +
  geom_point(aes(size = value), fill='#cca8e9',shape = 21, color = "black") + 
  labs(x = "", y = "", fill = "") +
  scale_size_continuous(breaks = c(50, 10000, 50000, 150000, 300000), range = c(1, 5))+
  expand_limits(x = x_range, y = y_range) +
  theme_classic()+
  theme(
    axis.text = element_text(size = 30, family = "Arial", color = "black"),    
    legend.text = element_text(size = 25, family = "Arial"),   
    legend.title = element_text(size = 25, family = "Arial"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),      
    axis.line = element_line(linewidth = 1.2, color = "black")        
  )
scatter_plot
# Save the heatmap with circles
ggsave('gatkSTR_heatmap_circles.png', plot = scatter_plot, dpi = 600,width = 8,height = 8)
############################################
library(ggplot2)

# Clear the environment
#rm(list = ls())

# Set the working directory

# Read the data from the file

# Create a heatmap with circles
scatter_plot <- ggplot(data_hip, aes(x = refrptime, y = hiprepeattime)) +
  geom_point(aes(size = value), fill='#f08a5d',shape = 21, color = "black") + 
  
  labs(x = "", y = "", fill = "") +
  expand_limits(x = x_range, y = y_range) +  
  scale_size_continuous(breaks = c(50, 10000, 50000, 150000, 300000), range = c(1, 5))+
  theme_classic()+
  theme(
    axis.text = element_text(size = 30, family = "Arial", color = "black"),      
    legend.text = element_text(size = 25, family = "Arial"),    
    legend.title = element_text(size = 25, family = "Arial"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),      
    axis.line = element_line(linewidth = 1.2, color = "black")        
  )
scatter_plot
# Save the heatmap with circles
ggsave('hipSTR_heatmap_circles.png', plot = scatter_plot, dpi = 600,width = 8,height = 8)
############################################################################################
scatter_plot <- ggplot(data_lob, aes(x = refrptime, y = lobrepeattime)) +
  geom_point(aes(size = value), fill='#95e1d3',shape = 21, color = "black") + 
  
  labs(x = "", y = "", fill = "") +
  expand_limits(x = x_range, y = y_range) +  
  scale_size_continuous(breaks = c(50, 10000, 50000, 150000, 300000), range = c(1, 5))+
  theme_classic()+
  theme(-
    axis.text = element_text(size = 30, family = "Arial", color = "black"),      
    legend.text = element_text(size = 25, family = "Arial"),    
    legend.title = element_text(size = 25, family = "Arial"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),      
    axis.line = element_line(linewidth = 1.2, color = "black")        
  )
scatter_plot
# Save the heatmap with circles
ggsave('lobSTR_heatmap_circles.png', plot = scatter_plot, dpi = 600,width = 8,height = 8)
