library(ggplot2)
library(data.table)
ns <- c(  5000, 10000, 500, 1000, 2500)
ns <- sort(ns)
hard_list <-   c(T,F)
pos_list <-  c(T,F)
use_oracle_sieve <- F
append <- "xgboost_2"
for(pos in pos_list){
  for(hard in hard_list) {
    try({
      sims_list <- lapply(ns, function(n) {
        load(paste0("mainSimResults/simsCATE", hard, pos,  "n", n, "_", append))
        simresults <- get(paste0("simresults", n))
        onestepbenchoracle <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbenchoracle")))
        onestepbench  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbench")))

        causalforestrisks <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_cf")))
        substrisks  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))



        lrnr_names <- simresults[[1]]$sieve$sieve_names #names(simresults[[1]]$CATEonestepbench) #simresults[[1]]$sieve[[1]]
        #lrnr_names <- unlist(lapply(lrnr_names, function(name) {
         # paste0(name ,"_", c( paste0("fourier_basis_", 0:4, "_plugin")))
        #}))
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
        dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier.basis_([0-9]+)")[,2])

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
         # tmp <- dt[, cv_sieve_risk := risks_oracle[which.min(cvrisksDR)], by = c("lrnr", "iter")]
        # tmp <- tmp[, oracle_sieve_risk := risks_oracle[which.min(risks_oracle)], by = c("lrnr", "iter")]
        # tmp <- tmp[!duplicated(paste0(degree, lrnr, iter )),]
        #
        #
        # tmp <- dt[, cv_sieve_risk := which.min(cvrisksDR), by = c("lrnr", "iter")]
        # tmp <- tmp[, oracle_sieve_risk := which.min(risks_oracle), by = c("lrnr", "iter")]
        # tmp <- tmp[!duplicated(paste0(degree, lrnr, iter )),]
        #


        ###### LATER


        #dt <- data.table(lrnr_full = lrnr_names, cvrisksDRoracle,cvrisksDR,  risks_oracle)
        # dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[,2])
        #dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0
        tmp <- data.table(risks_oracle = onestepbench, lrnr = names(onestepbench), lrnr_full = names(onestepbench), type = "DR", degree = "DR")
        dt <- rbind(dt, tmp, fill = T)

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

        unique(dt$lrnr)

        dt$type[!is.na(as.numeric(dt$degree))] <- "sieve"
        dt$type[is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]


        dt[, risks_best := mean(risks_oracle), by = c("lrnr", "type", "degree")]

        #dt[is.na(as.numeric(dt$degree)), risks_best := risks_oracle, by = c("lrnr", "type")]

        dt2 <- dt[,c("lrnr", "risks_best", "type", "degree"), with = F]
        dt2 <- unique(dt2)
        dt2$n <- n
        return(dt2)
      })


      dt <- rbindlist(sims_list)
     # dt <- dt[dt$type == "sieve",]
      dt <- dt[-grep("cv", dt$lrnr),]

      dt_tmp<-dt
      dt_tmp <- dt_tmp[grep("xgboost", dt$lrnr), ]
      dt_tmp$n <- as.factor(dt_tmp$n)
      max_depth <- stringr::str_match(dt_tmp$lrnr, "[0-9]+$")
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("xgboost (", "max_depth=", max_depth,")")

      dt_tmp_sieve <- dt_tmp[type == "sieve"]
      dt_tmp_sieve$degree <- as.numeric(dt_tmp_sieve$degree)
       plt <- ggplot(dt_tmp_sieve, aes(x = degree, y = risks_best, group = n, color = n, linetype = n)) + geom_line() +
        facet_wrap(~lrnr, scales = "free") + theme(axis.text.x = element_text(  vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = range(dt_tmp$risks_best))
      plt <- plt + geom_hline(data = dt_tmp[type != "sieve"], alpha = 0.4, aes(yintercept = risks_best,linetype= n, color = n))
      plt <-plt +  labs(group = "Sample Size (n)", color ="Sample Size (n)",  linetype = "Sample Size (n)",  "Order of trig. basis", y = "Mean-Squared-Error (MSE)" )


        ggsave(paste0("mainSimResults/SievePlotxgboost_oracleSieve", "pos=",pos, "hard=",hard, ".pdf"), width = 9, height = 6)



        dt_tmp<-dt[!(dt$lrnr %in% c("glm", "earth", "gam3", "gam4", "gam5")),]
        dt_tmp <- dt_tmp[-grep("xgboost", dt$lrnr), ]
        dt_tmp$n <- as.factor(dt_tmp$n)

        max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_xg$")[,2]
        max_depth[is.na(max_depth)] <- "cv"
        dt_tmp$lrnr <- paste0("ranger (", "max_depth=", max_depth,")")
        dt_tmp$lrnr <- as.factor(dt_tmp$lrnr)
        levels( dt_tmp$lrnr) <-  c(sort(unique(dt_tmp$lrnr[-grep("=[0-9][0-9]", dt_tmp$lrnr)])),
                     sort(unique(dt_tmp$lrnr[grep("=[0-9][0-9]", dt_tmp$lrnr)])))

        dt_tmp_sieve <- dt_tmp[type == "sieve"]
        dt_tmp_sieve$degree <- as.numeric(dt_tmp_sieve$degree)

        plt <- ggplot(dt_tmp_sieve, aes(x = degree, y = risks_best, group = n, color = n, linetype = n)) + geom_line() +
          facet_wrap(~lrnr, scales = "free") + theme(axis.text.x = element_text(  vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = range(dt_tmp$risks_best))
        plt <- plt + geom_hline(data = dt_tmp[type != "sieve"], alpha = 0.4, aes(yintercept = risks_best, linetype = n, color = n))
        plt <-plt +  labs(group = "Sample Size (n)", color ="Sample Size (n)",  linetype = "Sample Size (n)",  "Order of trig. basis", y = "Mean-Squared-Error (MSE)" )


        ggsave(paste0("mainSimResults/SievePlotRanger_oracleSieve", "pos=",pos, "hard=",hard, ".pdf"), width = 9, height = 6)

    })
  }
}
