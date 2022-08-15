#' confusion_matrix
#'
#' plot a confusion matrix based on two vectors of equal length
#'
#' @param vector_a first vector, e.g. true label
#' @param vector_b second vector, e.g. predicted label
#' @param xlab_text x label
#' @param ylab_text y label
#' @param plot_absolute if FALSE, show the relative fraction of values for each cell
#' @param plot_numbers show count in fields
#' @param text_size size of text
#' @return Returns a confusion matrix as heatmap
#' @export
#' @examples
#' true_label = c('A','A','A','B','B','B')
#' predicted_label = c('A','B','B','B','B','B')
#' #percentage of the total observations that are in the respective field
#'jj_plot_confusion_matrix(vector_a = true_label, vector_b = predicted_label,
#'                         xlab_text = 'true value', ylab_text = 'predicted value', 
#'                         plot_numbers = T, text_size = 10)
#'# absolute
#' jj_plot_confusion_matrix(vector_a = true_label, vector_b = predicted_label,
#'                         xlab_text = 'true value', ylab_text = 'predicted value', 
#'                         plot_numbers = T, text_size = 10, plot_absolute = T)

jj_plot_confusion_matrix = function(vector_a, vector_b, 
                                 xlab_text='vector_a', ylab_text='vector_b',
                                 plot_absolute = F, plot_numbers = F, text_size=2){
  library(cowplot)
  predictions <- table(vector_a, vector_b)
  if(!plot_absolute){
    # percentage of vector a (i.e. which fraction of total counts are in the respective field)
    predictions <- round(100*predictions/sum(predictions))  
  }
  predictions <- as.data.frame(predictions)
  colnames(predictions)= c('Var1', 'Var2', 'Freq')
  predictions$Freq_plot = ifelse(predictions$Freq ==0, NA, predictions$Freq)
  p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) +
    geom_tile() + 
    xlab(xlab_text) +
    ylab(ylab_text) + 
    theme_cowplot() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  if(plot_absolute){
    p1 = p1 + scale_fill_gradient(name = "Value count", low = "#ffffc8", high = "#7d0025") 
  }else{
    #Percentage of predicted values (vector_b)
    p1 = p1 + scale_fill_gradient(name = "Percentage", low = "#ffffc8", high = "#7d0025") 
  }
  if(plot_numbers){
    #label gives warning since NA are introduced in step above -> ignore
    p1 = p1 + geom_text(aes(label = Freq_plot), size=text_size, na.rm=TRUE)
  }
  return(p1)
}
