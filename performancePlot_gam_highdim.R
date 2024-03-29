library(ggplot2)
library(data.table)
ns <- c(  5000, 10000, 500, 1000, 2500)
ns <- sort(ns)
hard_list <-   c(T,F)
pos_list <-  c(T,F)
use_oracle_sieve <- F
for(pos in pos_list){
  for(hard in hard_list) {
    try({
    sims_list <- lapply(ns, function(n) {
      load(paste0("mainSimResults/simsCATE", hard, pos,  "n", n, "_gam_highDim"))
      simresults <- get(paste0("simresults", n))
      onestepbenchoracle <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbenchoracle")))
      onestepbench  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbench")))

      causalforestrisks <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_cf")))
      substrisks  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))



      lrnr_names <- names(simresults[[1]]$CATEonestepbench) #simresults[[1]]$sieve[[1]]
      lrnr_names <- unlist(lapply(lrnr_names, function(name) {
        paste0(name , c( "_no_sieve.plugin", paste0("_fourier_basis_", 1:4, "_plugin")))
      }))
      iter <- rep(1:length(simresults), each = length(lrnr_names))
      cvrisksDRoracle <- unlist( lapply(simresults, function(item) {
        item$sieve$cvrisksDRoracle
      }))

      cvrisksDR <- unlist(lapply(seq_along(simresults), function(index) {
        item <- simresults[[index]]
        as.vector(item$sieve$cvrisksDR)

      }))


      risks_oracle <- unlist( lapply(simresults, function(item) {
        item$sieve$risks_oracle
      }))

      dt <- data.table(iter, lrnr_full = lrnr_names, cvrisksDR, cvrisksDRoracle, risks_oracle)
      dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[,2])
      dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0

      dt$lrnr[ grep("gam3", dt$lrnr_full)] <- "gam3"
      dt$lrnr[ grep("gam4", dt$lrnr_full)] <- "gam4"
      dt$lrnr[ grep("gam5", dt$lrnr_full)] <- "gam5"
      dt$lrnr[ grep("glm", dt$lrnr_full)] <- "glm"
      dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
      dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
      dt$lrnr[ grep("rpart", dt$lrnr_full)] <- "rpart"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_7", dt$lrnr_full)] <- "ranger_7"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_13", dt$lrnr_full)] <- "ranger_13"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_10", dt$lrnr_full)] <- "ranger_10"

      ## dt$lrnr[ grep("xgboost_20_1_7", dt$lrnr_full)] <- "xgboost_7"
      #  dt$lrnr[ grep("xgboost_20_1_5", dt$lrnr_full)] <- "xgboost_5"
      # dt$lrnr[ grep("xgboost_20_1_3", dt$lrnr_full)] <- "xgboost_3"
      dt$lrnr <- gsub("Lrnr_", "", dt$lrnr)
      dt$lrnr <- gsub("_fourier_basis.+", "", dt$lrnr)
      dt$lrnr <- gsub(".fourier_basis.+", "", dt$lrnr)
      dt$lrnr <- gsub("_no_sieve.+", "", dt$lrnr)
      dt$lrnr <- gsub(".no_sieve.+", "", dt$lrnr)
      dt$lrnr <- gsub("_fourier.basis.+", "", dt$lrnr)
      dt$lrnr <- gsub(".fourier.basis.+", "", dt$lrnr)


      dt$type[!is.na(as.numeric(dt$degree))] <- "sieve"
      dt$type[is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

      dt <- dt[dt$degree > 0]

      tmp <- dt[, cv_sieve_risk := risks_oracle[which.min(cvrisksDR)], by = c("lrnr", "iter")]
      tmp <- tmp[, oracle_sieve_risk := risks_oracle[which.min(risks_oracle)], by = c("lrnr", "iter")]
      tmp <- tmp[!duplicated(paste0(degree, lrnr, iter )),]


       tmp <- dt[, cv_sieve_risk := which.min(cvrisksDR), by = c("lrnr", "iter")]
       tmp <- tmp[, oracle_sieve_risk := which.min(risks_oracle), by = c("lrnr", "iter")]
       tmp <- tmp[!duplicated(paste0(degree, lrnr, iter )),]



      ###### LATER


      #dt <- data.table(lrnr_full = lrnr_names, cvrisksDRoracle,cvrisksDR,  risks_oracle)
     # dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[,2])
      #dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0
      tmp <- data.table(risks_oracle = onestepbench, lrnr_full = names(onestepbench), degree = "DR")
      dt <- rbind(dt, tmp, fill = T)
      tmp <- data.table(risks_oracle = onestepbenchoracle, lrnr_full = names(onestepbench), degree = "DRoracle")
      dt <- rbind(dt, tmp, fill = T)
      tmp <- data.table(risks_oracle = causalforestrisks, lrnr_full = "causalforest", lrnr = "causalforest", degree = "causalforest")
      dt <- rbind(dt, tmp, fill = T)
      tmp <- data.table(risks_oracle = substrisks, lrnr_full = names(onestepbench), degree = "subst")
      dt <- rbind(dt, tmp, fill = T)
      dt$lrnr <- dt$lrnr_full
      dt$lrnr[ grep("gam3", dt$lrnr_full)] <- "gam3"
      dt$lrnr[ grep("gam4", dt$lrnr_full)] <- "gam4"
      dt$lrnr[ grep("gam5", dt$lrnr_full)] <- "gam5"
      dt$lrnr[ grep("glm", dt$lrnr_full)] <- "glm"
      dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
      dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
      dt$lrnr[ grep("rpart", dt$lrnr_full)] <- "rpart"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_7", dt$lrnr_full)] <- "ranger_7"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_13", dt$lrnr_full)] <- "ranger_13"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_10", dt$lrnr_full)] <- "ranger_10"

      ## dt$lrnr[ grep("xgboost_20_1_7", dt$lrnr_full)] <- "xgboost_7"
      #  dt$lrnr[ grep("xgboost_20_1_5", dt$lrnr_full)] <- "xgboost_5"
      # dt$lrnr[ grep("xgboost_20_1_3", dt$lrnr_full)] <- "xgboost_3"

      dt$lrnr <- gsub("Lrnr_", "", dt$lrnr)
      dt$lrnr <- gsub("_fourier_basis.+", "", dt$lrnr)
      dt$lrnr <- gsub(".fourier_basis.+", "", dt$lrnr)
      dt$lrnr <- gsub("_no_sieve.+", "", dt$lrnr)
      dt$lrnr <- gsub(".no_sieve.+", "", dt$lrnr)
      dt$lrnr <- gsub("_fourier.basis.+", "", dt$lrnr)
      dt$lrnr <- gsub(".fourier.basis.+", "", dt$lrnr)

      dt$type[!is.na(as.numeric(dt$degree))] <- "sieve"
      dt$type[is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

      print(unique(dt$lrnr_full))

      if(!use_oracle_sieve){
        dt[!is.na(as.numeric(dt$degree)), risks_oracle := risks_oracle[which.min(cvrisksDR)], by = c("lrnr", "type", "iter")]

      } else {
        dt[!is.na(as.numeric(dt$degree)), risks_oracle := min(risks_oracle), by = c("iter", "lrnr", "type")]
      }
      dt[, risks_best := mean(risks_oracle), by = c("lrnr", "type")]

      #dt[is.na(as.numeric(dt$degree)), risks_best := risks_oracle, by = c("lrnr", "type")]

      dt2 <- dt[,c("lrnr", "risks_best", "type"), with = F]
      dt2 <- unique(dt2)
      dt2$n <- n
      return(dt2)
    })


    dt <- rbindlist(sims_list)




    dt <- rbindlist(sims_list)
    dt <- dt[dt$lrnr!="causalforest",]
    s <- stringr::str_match(dt$lrnr[grep("gam",dt$lrnr)], "s(.+)_")[,2]
    dt$lrnr[grep("gam",dt$lrnr)] <- paste0("gam (s=", s, ")")
    dt$lrnr[grep("gam",dt$lrnr)][is.na(s)] <- paste0("gam (s=cv)")
    dt_tmp<-dt



    plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) + geom_line() +
      facet_wrap(~lrnr, scales = "free") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = c(min(1e-1, min(dt_tmp$risks_best)), max(dt_tmp$risks_best)))  +  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000))
    plt <- plt + xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") + theme_bw() + labs(color = "Method", group = "Method", linetype = "Method")
    ggsave(paste0("mainSimResults/performancePlot_GAM_CATEhighDim_", "pos=",pos, "hard=",hard, ".pdf"), width = 8, height = 7)



})
  }
}
