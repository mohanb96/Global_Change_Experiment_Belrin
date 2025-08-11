
#Estimating mean and its 95%confidence interval
BootStrap_mean = function(response, data=df, target = treatment, n_perm = n_iter){
  summary = list()
  
  for(treatment in target){
    bs = numeric(0)
    if(treatment=="1") population = data[data$remark%in%stressors, response]
    if(treatment!="1") population = data[data$remark==treatment, response]
    size = length(population)-sum(is.na(population))
    
    for(id in c(1:n_perm)){
      k = mean(sample(population, size, replace = T), na.rm = TRUE)
      bs = append(bs, k)
    }
    summary[[treatment]] = c(quantile(bs, .025,na.rm = TRUE), mean(bs,na.rm = TRUE), quantile(bs, .975,na.rm = TRUE))
    names(summary[[treatment]]) = c("2.5%", "mean", "97.5%")
  }
  summary = t(data.frame(summary))
  summary = data.frame("target" = target, summary); row.names(summary) = c()
  return(summary)
}

BootStrap_ES_rep = function(response, data=df, target = treatment, n_perm = n_iter){
  resampled = list()
  
  population_CT = data[data$remark=="CT", response]
  
  for(treatment in target){
    bs = numeric(0)
    if(treatment=="1") population_TR = data[data$remark%in%stressors, response]
    if(treatment!="1") population_TR = data[data$remark==treatment, response]
    size_CT = length(population_CT)-sum(is.na(population_CT))
    size_TR = length(population_TR)-sum(is.na(population_TR))
    
    for(id in c(1:n_perm)){
      k_CT = mean(sample(population_CT, size_CT, replace = T), na.rm = TRUE)
      k_TR = mean(sample(population_TR, size_TR, replace = T), na.rm = TRUE)
      bs = append(bs, k_TR - k_CT)
    }
    resampled[[treatment]] = bs
  }
  resampled[["CT"]] = rep(0, n_perm)
  return(resampled)
}

BootStrap_ES_summary = function(data){
  summary = list()
  p = 0
  summary[["CT"]] = c(0,0,0,1)
  target = names(data)
  
  for(treatment in target[-1]){
    bs = data[[treatment]]
    p = length(which(bs>0))/length(bs)
    p = min(p, 1-p)
    summary[[treatment]] = c(quantile(bs, .025,na.rm = TRUE), mean(bs,na.rm = TRUE), quantile(bs, .975,na.rm = TRUE), p)
  }
  summary = t(data.frame(summary))
  colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
  summary = data.frame(target, summary); row.names(summary) = c()
  
  return(summary)
}

Null_distribution_rep = function(response, data=df, n_perm=n_iter){
  
  output = list()
  for(Lv in levels){
    
    resampled = list()
    
    # Checking which stressor combinations were jointly tested
    if(Lv=="1") combination = data[data$remark%in%stressors,stressors]
    if(Lv!="1") combination = data[data["remark"]==Lv,stressors]
    Level = sum(combination[1,]) 
    
    
    # Null distributions can be taken based on three different assumptions
    for(type in c("Additive", "Multiplicative", "Dominative")){
      
      population_CT = df[df$remark=="CT", response]
      size_CT = length(population_CT)-sum(is.na(population_CT)) ##subtract NA value
      
      # For each combination, bootstrap resampling is conducted
      for(j in c(1:nrow(combination))){
        bs = numeric(0)
        selected_stressors = stressors[which(combination[j,]==1)]
        sub_n_perm = ceiling(n_perm/nrow(combination)) #*5 increase the permutation number for sub-sampling
        
        # bootstrap resampling
        for(id in c(1:sub_n_perm)){
          each_effect = numeric(0)
          k_CT = mean(sample(population_CT, size_CT, replace = T),na.rm = TRUE) #ignore NA
          
          for(treatment in selected_stressors){
            population_TR = df[df$remark==treatment, response]
            size_TR = length(population_TR)-sum(is.na(population_TR))
            k_TR = mean(sample(population_TR, size_TR, replace = T),na.rm = TRUE)#ignore NA
            
            # ES estimate depending on the type of null hypotheses
            if(type=="Additive")       each_effect = append(each_effect, (k_TR - k_CT))
            if(type=="Multiplicative") each_effect = append(each_effect, (k_TR - k_CT)/k_CT)
            if(type=="Dominative")      each_effect = append(each_effect, (k_TR - k_CT))
          }
          
          # Calculating an expected ES after collecting the ESs of all relevant single stressors
          if(type=="Additive")       joint_effect = sum(each_effect)
          if(type=="Multiplicative"){
            z = 1
            for(m in c(1:Level)) z = z * (1 + each_effect[m])
            joint_effect = (z - 1)*k_CT
          }
          if(type=="Dominative")      joint_effect = each_effect[which(max(abs(each_effect))==abs(each_effect))]
          
          bs = append(bs, joint_effect)
        }
        resampled[[type]][[j]] = bs
      }
      
    }
    output[[Lv]] = resampled
  }  
  return(output)
}



Null_distribution_rep_transform = function(data){
  output = list()
  for(Lv in levels){
    for(type in c("Additive", "Multiplicative", "Dominative")){
      output[[Lv]][[type]] = sample(unlist(data[[Lv]][[type]]), n_iter, replace=F)
    }
  }
  return(output)
}

NHST_summary = function(null_data, Actual_data){
  output = list()
  for(Lv in levels){
    summary = list()
    summary[["Actual"]] = c(quantile(Actual_data[[Lv]], .025,na.rm = TRUE), mean(Actual_data[[Lv]],na.rm = TRUE), quantile(Actual_data[[Lv]], .975,na.rm = TRUE), 1)
    p = 0
    assumptions = c("Additive", "Multiplicative", "Dominative")
    
    for(i_assumption in assumptions){
      bs   = (Actual_data[[Lv]] - null_data[[Lv]][[i_assumption]])
      p = length(which(bs>0))/length(bs)
      p = min(p, 1-p)
      summary[[i_assumption]] = c(quantile(null_data[[Lv]][[i_assumption]], .025,na.rm = TRUE), mean(null_data[[Lv]][[i_assumption]],na.rm = TRUE), quantile(null_data[[Lv]][[i_assumption]], .975,na.rm = TRUE), p)
    }
    summary = t(data.frame(summary))
    colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
    summary = data.frame(ES = c("Actual", "Additive","Multiplicative","Dominative"), summary); row.names(summary) = c()
    
    output[[Lv]] = summary
  }
  
  return(output)
}

NHST_summary_transform = function(data){
  output = list()
  for(i in 1:4){
    summary = rbind(data[["1"]][i, 2:4], data[["2"]][i, 2:4], data[["4"]][i, 2:4],
                    data[["6"]][i, 2:4], data[["8"]][i, 2:4], data[["10"]][i, 2:4])
    summary = cbind(levels, summary)
    colnames(summary) = c("Lv", "Low", "Mean", "High")
    output[[c("Actual", "Additive", "Multiplicative", "Dominative")[i]]] = summary
  }
  return(output)
}

Expected_ES_for_each = function(data){
  output = numeric(0)
  for(type in c("Additive", "Multiplicative", "Dominative")){
    tmp = numeric(0)
    for(Lv in levels){
      n_len = length(data[[Lv]][[type]])
      for(i in 1:n_len){
        tmp = append(tmp, mean(data[[Lv]][[type]][[i]]))
      }
    }
    output = cbind(output,tmp)
  }
  colnames(output)= c("E1", "E2", "E3")
  return(output)
}


postResample <- function(pred, obs)
{
  isNA <- is.na(pred)
  pred <- pred[!isNA]
  obs <- obs[!isNA]
  if (!is.factor(obs) && is.numeric(obs))
  {
    if(length(obs) + length(pred) == 0)
    {
      out <- rep(NA, 3)
    } else {
      if(length(unique(pred)) < 2 || length(unique(obs)) < 2)
      {
        resamplCor <- NA
      } else {
        resamplCor <- try(cor(pred, obs, use = "pairwise.complete.obs"), silent = TRUE)
        if (inherits(resamplCor, "try-error")) resamplCor <- NA
      }
      mse <- mean((pred - obs)^2)
      mae <- mean(abs(pred - obs))
      out <- c(sqrt(mse), resamplCor^2, mae)
    }
    names(out) <- c("RMSE", "Rsquared", "MAE")
  } else {
    if(length(obs) + length(pred) == 0)
    {
      out <- rep(NA, 2)
    } else {
      pred <- factor(pred, levels = levels(obs))
      requireNamespaceQuietStop("e1071")
      out <- unlist(e1071::classAgreement(table(obs, pred)))[c("diag", "kappa")]
    }
    names(out) <- c("Accuracy", "Kappa")
  }
  if(any(is.nan(out))) out[is.nan(out)] <- NA
  out
}

### Specific null model prediction###
NullModel = function(response, data = df, selected_factors=vector(), n_perm = 100){
  
  output = list()
  population_CT= data[data$remark == "CT", response]
  size_CT=length(population_CT)
  CT = mean(population_CT)
  
  for (type in c("additive","multiplicative", "dominative")) {
    bs = numeric(0)
    for (id in c(1:n_perm)) {
      each_effect = numeric(0)
      k_CT = mean(sample(population_CT, size_CT, replace = T))
      for (treatment in selected_factors) {
        population_TR = data[data$remark == treatment, response]
        size_TR = length(population_TR)
        k_TR = mean(sample(population_TR, size_TR, replace = T))
        
        # ES estimate depending on the type of null hypothesis
        if(type == "additive")        each_effect = append(each_effect,(k_TR - k_CT))
        if(type == "multiplicative")  each_effect = append(each_effect, (k_TR-k_CT)/k_CT)
        if(type == "dominative")      each_effect = append(each_effect,(k_TR - k_CT))
      }
      
      if(type == "additive") {
        joint_effect = sum(each_effect)
        pre_response = joint_effect+CT
      }
      
      if(type=="multiplicative"){
        z = 1
        for(m in c(1:length(selected_factors))) {
          z = z * (1 + each_effect[m])
          joint_effect = (z - 1)*k_CT
        }
        pre_response =joint_effect+CT
      }
      
      if(type=="dominative")  {
        joint_effect = each_effect[which(max(abs(each_effect))==abs(each_effect))]
        pre_response =joint_effect+CT
      }
      
      bs = append(bs, pre_response)
    }
    output[[type]] = bs
  }
  return(output)
}

#Null_modle_i_Treatment = NullModel("Tp_cover_mar22", selected_factors = c("Warming","N_Addition"),n_perm = n_perm)

NullModel_summary = function(null_data, actual_data = vector()){
  output = list()
  
  if(length(actual_data)<=2){
    output[["actual"]] = c(mean(actual_data),mean(actual_data),mean(actual_data),1)
    for (type in c("additive","multiplicative", "dominative")) {  
      bs = rep(mean(actual_data),length(null_data[[type]])) - null_data[[type]]
      p = length(which(bs>0))/length(bs)
      p = min(p, 1-p)
      output[[type]] = c(quantile(null_data[[type]], .025), mean(null_data[[type]]), quantile(null_data[[type]], .975), p)
    }
  }
  
  if(length(actual_data)>=3){
    #re-sampling of actual data
    size_actul=length(actual_data)
    bs_actual = numeric(0)
    for (i in c(1:length(null_data[[1]]))) {
      k_actual = mean(sample(actual_data, size_actul, replace = T))
      bs_actual = append(bs_actual, k_actual)
    }
    
    output[["actual"]] = c(quantile(bs_actual, .025), mean(bs_actual), quantile(bs_actual, .975), 1)
    
    for (type in c("additive","multiplicative", "dominative")) { 
      bs = bs_actual - null_data[[type]]
      p = length(which(bs>0))/length(bs)
      p = min(p, 1-p)
      output[[type]] = c(quantile(null_data[[type]], .025), mean(null_data[[type]]), quantile(null_data[[type]], .975), p)
    }
  }
  
  output = t(data.frame(output))
  colnames(output) = c("X2.5%", "mean", "X97.5%", "p_value")
  output = data.frame(ES = c("actual", "additive","multiplicative","dominative"), output)
  row.names(output) = c()
  
  return(output)
}

###This function is used for microbial dataset transformation.
micro_data_transform <- function(data, phylums, meta_data = Null){
  
  phylum_list <- list()
  
  for (phylum_i in phylums) {
    
    phylum_list[[phylum_i]] <- c()
    data_phylum_i <- data[data[,"phylum"]==phylum_i,]
    
    if (nrow(data_phylum_i) < 154){
      
      ##Add NA value to the missing plots
      existing_ids <- data_phylum_i$Plot_ID
      missing_ids <- setdiff(1:154, existing_ids)
      missing_df <- data.frame(Plot_ID = missing_ids)
      missing_df[["relative_abundance"]] <- NA
      other_cols <- setdiff(names(data_phylum_i), c("Plot_ID","relative_abundance"))
      missing_df <- missing_df%>%
        mutate(phylum = rep(phylum_i,nrow(missing_df)))%>%
        mutate(taxonomic_group = rep(data_phylum_i[1,"taxonomic_group"],nrow(missing_df)))
      data_phylum_i <- rbind(data_phylum_i,missing_df)
      data_phylum_i <- arrange(data_phylum_i,by = Plot_ID)
    }
    
    for (id in 1:154) {
      phylum_list[[phylum_i]][id] <- data_phylum_i[data_phylum_i[,"Plot_ID"]==id,"relative_abundance"]
    }
  }
  
  Meta_data <- meta_data
  for (phylum_i in phylums) {
    Meta_data[[phylum_i]] <- phylum_list[[phylum_i]]
  }
  
  return(Meta_data)
}
