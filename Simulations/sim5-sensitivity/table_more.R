library(here)
library(huxtable)
setwd(here('sim5-sensitivity'))

## CATE
cases <- c('binaryY-saturated','continuousY-binaryX','continuousY-continuousX')

dgps <- c('fd_admg1','fd_admg2','nonfd_admg1','nonfd_admg2')

dat.table1 <- c(500,1000,2000,4000,10000)
# dat.table1 <- data.frame(N=c(500,1000,2000,4000),V1=rep(NA,4),V2=rep(NA,4),V3=rep(NA,4),V4=rep(NA,4))

for (i in 1:length(cases)) {
  
  dat <- rep(1,5)
  
  for (j in 1:length(dgps)) {
    load(paste0(cases[i],'/',dgps[j],'/result_AIPW.Rdata'))
    
    dat <- cbind(dat,cate.wald.rate_matrix[,2])
    
  }
  
  dat <- dat[,-1]
  
  dat.table1 <- cbind(dat.table1,dat)
}

# dat.table1 <- cbind(dat.table1,data.frame(V1=rep(NA,5),V2=rep(NA,5),V3=rep(NA,5),V4=rep(NA,5)))

colnames(dat.table1) <- c('N',paste0('V',1:12))

cases <- c('binaryY-saturated','continuousY-binaryX','continuousY-continuousX')

dgps <- c('fd_admg1','fd_admg2','nonfd_admg1','nonfd_admg2')


dat.table2 <- c(500,1000,2000,4000,10000)

for (i in 1:length(cases)) {
  
  dat <- rep(1,5)
  
  for (j in 1:length(dgps)) {
    load(paste0(cases[i],'/',dgps[j],'/result_AIPW.Rdata'))
    
    dat <- cbind(dat,dual.wald.rate_matrix[,2])
    
  }
  
  dat <- dat[,-1]
  
  dat.table2 <- cbind(dat.table2,dat)
}

colnames(dat.table2) <- c('N',paste0('V',1:12))

dat.table3 <- c(500,1000,2000,4000,10000)
for (i in 1:length(cases)) {
  
  dat <- rep(1,5)
  
  for (j in 1:length(dgps)) {
    load(paste0(cases[i],'/',dgps[j],'/result_AIPW.Rdata'))
    
    dat <- cbind(dat,primal.wald.rate_matrix[,2])
    
  }
  
  dat <- dat[,-1]
  
  dat.table3 <- cbind(dat.table3,dat)
}

colnames(dat.table3) <- c('N',paste0('V',1:12))

dat.table <- rbind(dat.table1,dat.table2,dat.table3)

# make table
library(huxtable)
library(dplyr)
library(stringr)
options(scipen = 999)

dat.table <- as.matrix(dat.table)

table1 <- as_hux(dat.table[,-1]) %>% 
  insert_row("All binary","","","","Y continuous","","","","Y,X continuous","","","", after = 0) %>% 
  insert_row("Type I error","","Power","","Type I error","","Power","","Type I error","","Power","", after = 1) %>% 
  insert_row("ADMG1","ADMG2","ADMG3","ADMG4","ADMG1","ADMG2","ADMG3","ADMG4","ADMG1","ADMG2","ADMG3","ADMG4",after=2) %>%
  insert_row("CCM test","","","","","","","","","","","",after=3) %>%
  insert_row("Dual test","","","","","","","","","","","",after=9) %>%
  insert_row("Primal test","","","","","","","","","","","",after=15) %>%
  merge_cells(1,1:4) %>% 
  merge_cells(1,5:8) %>%
  merge_cells(1,9:12) %>%
  merge_cells(2, 1:2) %>% 
  merge_cells(2, 3:4) %>% 
  merge_cells(2, 5:6) %>% 
  merge_cells(2, 9:10) %>% 
  merge_cells(2, 11:12) %>% 
  # merge_cells(2, 13:14) %>% 
  insert_column(c("","N","","CATE-based test","500","1000","2000","4000","10000","Dual test","500","1000","2000","4000","10000","Primal test","500","1000","2000","4000","10000"), after = 0) %>%
  # set_number_format(6, 1, function(x) format(x, scientific = FALSE)) %>% 
  set_number_format(everywhere,everywhere, "%s") %>% 
  merge_cells(4, 1:3) %>%
  merge_cells(10, 1:3) %>%
  merge_cells(16, 1:3) %>%
  set_align(col=1, everywhere, "left") %>%
  set_align(col=2:ncol(.),everywhere,"center") %>%
  set_all_padding(-2)  %>% 
  set_bold(1, everywhere) %>%
  # set_bold(c(3,9,15), everywhere) %>%
  set_italic(2,everywhere) %>%
  set_bottom_border(row = 1, col =2:ncol(.)) %>% 
  # set_bottom_border(row = 3, col =2:ncol(.)) %>% 
  set_bottom_border(row = 2, col =2:ncol(.)) %>% 
  set_bottom_border(row = 3, col =2:ncol(.)) %>% 
  # set_right_border(4:nrow(.), 4, brdr(0.4,"solid")) %>%
  set_right_border(4:nrow(.), 5, brdr(0.4,"double")) %>%
  set_right_border(4:nrow(.), 9, brdr(0.4,"double")) %>%
  # set_right_border(4:nrow(.), 13, brdr(0.4,"double")) %>%
  # set_right_border(4:nrow(.), 10, brdr(0.4,"solid")) %>%
  # set_right_border(4:nrow(.), 13, brdr(0.4,"double")) %>%
  # set_right_border(4:nrow(.), 16, brdr(0.4,"solid")) %>%
  set_top_border(row=1,col=everywhere,brdr(1, "solid")) %>% set_bottom_border(row = nrow(.),col = everywhere,brdr(1, "solid")) %>%
  set_font_size(6) %>% set_caption("Evaluation of CATE-based, dual, and primal tests on type I error and power across various data settings.")

table1

quick_latex(table1, file = "table_more.tex")
