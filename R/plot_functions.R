#----------------------------------------------------------------------------------------------------
#' Plot ROC Curve
#'
#' Plot a ROC curve for the given statistics dataframe.
#'
#' @param stat.dataframe A data frame containing statistics for a particular model or group of
#' models. It must contain columns for "specificity", "sensitivity", and "Model.Name", else an error
#' will be returned
#' @param printAUC A Boolean variable indicating whether or not the function should print the
#' area under the curve on the plot (default = FALSE)
#'
#' @return A plot of ROC curves for all models in the statistics dataframe
#'
#' @export

plot.roc.curve <- function(stat.dataframe, printAUC=FALSE) {
    
    if (printAUC == TRUE){
        specif.sense.df <- stat.dataframe[,c("specificity","sensitivity")]
        ordered.specif.sense.df <- specif.sense.df[with(specif.sense.df, order(specificity)), ]
        ss.x <- ordered.specif.sense.df$specificity
        ss.y <- ordered.specif.sense.df$sensitivity
        AUC <- caTools::trapz(ss.x,ss.y)
    }
    
    p <- ggplot2::ggplot(stat.dataframe) +
        ggplot2::geom_line(ggplot2::aes(x=specificity, y=sensitivity, color=Model.Name)) +
        ggplot2::geom_abline(intercept = 1, slope = 1, color='lightgrey') +
        ggplot2::scale_x_reverse(lim=c(1,0)) +
        ggplot2::scale_y_continuous(lim=c(0,1)) +
        ggplot2::coord_fixed(ratio=1) +
        ggplot2::labs(x="Specificity", y="Sensitivity") +
        ggplot2::theme_minimal()
    if (printAUC == TRUE) {
        p <- p + ggplot2::annotate("text", x=0.25, y=0.2,
                                   label=paste('AUC = ', format(AUC,digits=5), sep=''))
    }
    p
} #plot.roc.curve
#----------------------------------------------------------------------------------------------------
#' Plot Precision-Recall Curve
#'
#' Plot a Precision-Recall curve for the given statistics dataframe
#'
#' @param stat.dataframe A data frame containing statistics for a particular model or group of
#' models. It must contain columns for "sensitivity", "ppv", and "Model.Name"
#'
#' @return A plot of Precision-Recall curves for all models in the statistics dataframe
#'
#' @export

plot.precrecall.curve <- function(stat.dataframe) {
    ggplot2::ggplot(stat.dataframe) + 
        ggplot2::geom_line(ggplot2::aes(x=sensitivity, y=ppv, color=Model.Name)) + 
            ggplot2::scale_x_continuous(lim=c(0,1)) + ggplot2::scale_y_continuous(lim=c(0,1)) +
            ggplot2::coord_fixed(ratio=1) +
            ggplot2::labs(x="Precision", y="Recall") +
            ggplot2::theme_minimal()
    
} #plot.precrecall.curve
#----------------------------------------------------------------------------------------------------
#' Plot a Matthews Correlation Coefficient Curve
#'
#' Plot a Matthews Correlation Coefficient Curve for the given statistics dataframe
#'
#' @param stat.dataframe A data frame containing statistics for a particular model or group of
#' models. It must contain columns for "threshold", "MattCC", and "Model.Name"
#'
#' @return A plot of Matthews Correlations Coefficient curves for all models in the statistics
#' data frame
#'
#' @export

plot.mattcc.curve <- function(stat.dataframe) {
    ggplot2::ggplot(stat.dataframe) +
        ggplot2::geom_line(ggplot2::aes(x=threshold, y=MattCC, color=Model.Name)) +
            ggplot2::scale_x_continuous(lim=c(0,1)) + ggplot2::scale_y_continuous(lim=c(0,1)) +
            ggplot2::coord_fixed(ratio=1) +
            ggplot2::labs(x="Threshold", y="Matthews correlation coefficient") +
            ggplot2::theme_minimal()

} #plot.mattcc.curve
#----------------------------------------------------------------------------------------------------
#' Plot All Curves
#'
#' Plot all parameters from the statistics dataframe
#'
#' @param stat.dataframe A data frame containing statistics for a particular model or group of
#' models. It must contain the column "threshold".
#'
#' @return A plot of all curves, where each curve corresponds to a different statistic. Note that
#' for this plot, the statistics are not separated by model type.
#'
#' @export

plot.all.curve <- function(stat.dataframe) {
    stat.dataframe.tidy <- tidyr::gather(stat.dataframe,
                                         'Statistic',
                                         'Value',
                                         -threshold,
                                         -Model.Name)    
    ggplot2::ggplot(stat.dataframe.tidy) +
        ggplot2::geom_line(ggplot2::aes(x=threshold, y=Value, color=Statistic)) +
            ggplot2::scale_x_continuous(lim=c(0,1)) + ggplot2::scale_y_continuous(lim=c(0,1)) +
            ggplot2::coord_fixed(ratio=1) +
            ggplot2::labs(x="Threshold") +
            ggplot2::theme_minimal()
}


