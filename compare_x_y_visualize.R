library(cowplot)
library(ggplot2)



keyx <- "PARfiltered.chrX"
keyy <- "PARfiltered.chrY"


  
  X_count <- read.csv(paste(keyx,"_changed.variety.tc_0815.csv", sep=""))
  Y_count <- read.csv(paste(keyy,"_changed.variety.tc_0815.csv",sep=""))
  X_count[is.na(X_count)] <- 0
  Y_count[is.na(Y_count)] <- 0
  rownames(X_count) <- X_count[, 1]  ## set rownames
  X_count <- X_count[, -1] 
  head(X_count)
  rownames(Y_count) <- Y_count[, 1]  ## set rownames
  Y_count <- Y_count[, -1] 
  head(Y_count)
  
  table_a <- as.matrix(X_count)
  table_b <- as.matrix(Y_count)
  singlebase_barplot_list <- list()
  new_table_a <- table_a
  new_table_b <- table_b
  rownames(table_a)
  
  current_rownames <- rownames(table_a)
  print(current_rownames)

 
  rownames(new_table_a) <- substr(current_rownames,1,nchar(current_rownames)-4)
  rownames(new_table_b) <- substr(current_rownames,1,nchar(current_rownames)-4)
  table_a <- rowsum(new_table_a, row.names(new_table_a))
  table_b <- rowsum(new_table_b, row.names(new_table_b))
  column_names <- colnames(table_a)
  
  

  XdifY <- data.frame(matrix(ncol = 10, nrow = 6))
  chrX <- data.frame(matrix(ncol = 10, nrow = 6))
  chrY <- data.frame(matrix(ncol = 10, nrow = 6))
  
  chrX_counts <- data.frame(matrix(ncol = 10, nrow = 6))
  chrY_counts <- data.frame(matrix(ncol = 10, nrow = 6))
  combined_counts <- list()
  combined <- list()
  for (i in 1:10) {
    column_name <- column_names[[i+2]]
    print(column_name)
    column_a <- table_a[, column_name]
    column_b <- table_b[, column_name]
    print(column_a)
    print(sum(column_a))

    
    combined_counts[[i]] <- data.frame(
      RowNames = rownames(table_a),
      chrX_counts = column_a,
      chrY_counts = column_b
      
    )
    
    combined[[i]] <- data.frame(
      RowNames = rownames(table_a),
      chrX = column_a/sum(column_a),
      chrY = column_b/sum(column_b)
    )
    print(combined[[i]])

    combined_long <- tidyr::pivot_longer(combined[[i]], -RowNames, names_to = "Table", values_to = "Values")
    print(combined_long)


    print(combined[[i]]$XoverY)
    combined[[i]]$XdifY <- combined[[i]]$chrX-combined[[i]]$chrY
    
   
    XdifY[,i] <- combined[[i]]$XdifY
    chrX[,i] <- combined[[i]]$chrX
    chrY[,i] <- combined[[i]]$chrY
    chrX_counts[,i] <- combined_counts[[i]]$chrX_counts
    chrY_counts[,i] <- combined_counts[[i]]$chrY_counts
  }
  

  
  XdifY_out <- paste(keyx,".six_mutation_xdify_0322.csv", sep="")
  colnames(XdifY) <- column_names[3:12]
  rownames(XdifY) <- rownames(combined[[1]])
  XdifY$p <- c(t.test(XdifY[1,],rep(0,10))$p.value,t.test(XdifY[2,],rep(0,10))$p.value,t.test(XdifY[3,],rep(0,10))$p.value,t.test(XdifY[4,],rep(0,10))$p.value,t.test(XdifY[5,],rep(0,10))$p.value,t.test(XdifY[6,],rep(0,10))$p.value)
  XdifY$corrected.p <- p.adjust(XdifY$p, method = "bonferroni")
  write.csv(XdifY, XdifY_out)
  
  chrX_out <- paste(keyx,".six_mutation_chrX.csv", sep="")
  chrY_out <- paste(keyx,".six_mutation_chrY.csv", sep="")
  
  write.csv(chrX, chrX_out)
  write.csv(chrY, chrY_out)


colnames(chrX) <- column_names[3:12]
rownames(chrX) <- rownames(combined[[1]])
colnames(chrY) <- column_names[3:12]
rownames(chrY) <- rownames(combined[[1]])

x_mutation <- t(chrX)
y_mutation <- t(chrY)

rownames(x_mutation) <- rep("chrX", 10)
rownames(y_mutation) <- rep("chrY", 10)


mutation <- rbind(x_mutation, y_mutation)
print(head(mutation))
reshape2::melt(mutation, variable.name = "mutation") -> mutation_wide
mutation_wide$Var2 <- gsub("(.)(.)$", "\\1>\\2", mutation_wide$Var2)

mutation_wide_out <- paste("X_Y",".six_mutation_wide_0322.csv", sep="")
write.csv(mutation_wide,mutation_wide_out, row.names = FALSE)
print(head(mutation_wide))

mutation_wide$Var1 <- factor(mutation_wide$Var1, levels = c("chrX", "chrY"))



## only plot chrx and chry


p_box_nopar <- ggplot(mutation_wide,aes(Var2,value)) +
  geom_boxplot(aes(fill=Var1),position = position_dodge(width = 0.8),size=0.3,outlier.shape = NA,width=0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8), size = 2, aes(fill = Var1),shape=21,alpha=0.3) + # Add this line for points
  theme_bw(base_size=14) +
  
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "grey30"),axis.text.x = element_text(angle =45, vjust = 0.5, hjust=0.6)) +
  theme(aspect.ratio=1,legend.position="top",
        panel.grid.major = element_line(color = "grey70", linetype = "longdash", size=0.1), 
        panel.grid.minor = element_blank(), ) +
  labs(x="\nType of substitutions",y="Frequency among total substitutions per branch\n",fill="Chromosome")

    

ggsave(plot =p_box_nopar,"single_mutation_dot_box_0322_with_p_nopar.pdf",width = 8,height = 5)
