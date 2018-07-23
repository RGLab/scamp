#copy of function "external_validation" in ClusterR package https://cran.r-project.org/web/packages/ClusterR/index.html
#with one modification: adjusted to provide methodName parameter
#this version is not affiliated with the ClusterR package or its auther/maintainer.
myExternalValidation = function(true_labels, clusters, methodName) {
    #modified from ClusterR package so it returns a column vector of clustering methods
    if (is.integer(true_labels)) true_labels = as.numeric(true_labels)
    if (is.integer(clusters)) clusters = as.numeric(clusters)
    if (!is.vector(true_labels) || !is.numeric(true_labels)) stop('true_labels should be a numeric vector')
    if (!is.vector(clusters) || !is.numeric(clusters)) stop('clusters should be a numeric vector')
    if (length(true_labels) != length(clusters)) stop('the length of the true_labels vector should equal the length of the clusters vector')

    tbl = table(clusters, true_labels)

    conv_df = as.data.frame.matrix(tbl)

                                        # Diagonal = rep(0, ncol(conv_df))
                                        #
                                        # for (i in 1:nrow(conv_df)) {
                                        #
                                        #   wh_idx = which.max(conv_df[i, ])
                                        #
                                        #   if (conv_df[i, wh_idx] > Diagonal[wh_idx]) {
                                        #
                                        #     Diagonal[wh_idx] = conv_df[i, wh_idx]
                                        #   }
                                        # }

    conv_df = as.data.frame.matrix(tbl)

    tp_plus_fp = sum(gmp::asNumeric(gmp::chooseZ(rowSums(conv_df), 2)))

    tp_plus_fn = sum(gmp::asNumeric(gmp::chooseZ(colSums(conv_df), 2)))

    tp = sum(gmp::asNumeric(gmp::chooseZ(as.vector(as.matrix(conv_df)), 2)))

    fp = tp_plus_fp - tp

    fn = tp_plus_fn - tp

    tn = gmp::asNumeric(gmp::chooseZ(sum(as.vector(as.matrix(conv_df))), 2)) - tp - fp - fn


    prod_comb = (tp_plus_fp * tp_plus_fn) / gmp::asNumeric(gmp::chooseZ(length(true_labels), 2))

    mean_comb = (tp_plus_fp + tp_plus_fn) / 2.0



    tmp_pur = apply(conv_df, 1, max)

    res_purity = sum(tmp_pur)/length(true_labels)


    tmp_entropy = sum(apply(conv_df, 2, function(x) entropyFormula(x)))

    res_entropy = -(1/(sum(tbl) * log2(length(unique(true_labels))))) * tmp_entropy


    mutual_information = 0.0

    joint_entropy = 0.0

    for (i in 1:nrow(conv_df)) {

        for (j in 1:ncol(conv_df)) {

            if (conv_df[i,j] > 0.0) {

                joint_entropy = joint_entropy + (-((conv_df[i,j] / sum(tbl)) * log2(conv_df[i,j] / sum(tbl))))

                mutual_information = mutual_information + ((conv_df[i,j] / sum(tbl)) * log2((sum(tbl) * conv_df[i,j]) / (sum(conv_df[i,]) * sum(conv_df[,j]))))
            }
        }
    }

    entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))

    entr_class = sum(apply(conv_df, 2, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))

    NMI = (mutual_information / ((entr_cluster + entr_class) / 2.0))

    VAR_INFO = (entr_cluster + entr_class) - 2.0 * mutual_information

    NVI = 1.0 - (mutual_information / joint_entropy)


    prec = tp / (tp + fp)
    rec = tp / (tp + fn)

    ovec <- c(round(res_purity, 3),
              round(res_entropy,3),
              round(NMI,3),
              round(VAR_INFO,3),
              round(NVI,3),
              round(tn / (tn + fp), 3),
              round(tp / (tp + fn), 3),
              round(prec,3),
              round(rec,3),
              round(2.0 * ((prec * rec) / (prec + rec)), 3),
              round((tp + tn) / (tp + fp + fn + tn), 3),
              round((tp - prod_comb) / (mean_comb - prod_comb), 3),
              round(tp / (tp + fp + fn), 3),
              round(sqrt((tp / ((tp + fp))) * (tp / (tp + fn))), 3))#,
              #round(2.0 * (fp + fn), 4))
    omat <- matrix(ovec,ncol=1)
    rownames(omat) <- c("Purity","Entropy","NormMutInfo","VarInfo","NormVarInfo","Specificity","Sensitivity",
                        "Precision","Recall","F-measure","Rand Index", "Adj Rand Index", "Jaccard Index",
                        "Fowlkes-Mallows Index")#, "Mirkin Metric")
    
    colnames(omat) <- methodName
    return(omat)
    
}


entropyFormula = function(x_vec) {
    #copied from ClusterR package. Used in modified external_validation call.
    vec = rep(NA, length(x_vec))
    for (i in 1:length(x <- vec)) {
        if (x_vec[i] == 0.0) {
            vec[i] = 0.0}
        else {
            vec[i] = ((x_vec[i]) * log2(x_vec[i]/sum(x_vec)))
        }
    }
    return(vec)
}

