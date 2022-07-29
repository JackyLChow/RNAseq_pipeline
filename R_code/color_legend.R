# color legend

lgd <- list()

for(m_s_ in srt){
  if(is.function(col_leg[[m_s_]])){
    lgd[[m_s_]] <- Legend(col_ = col_leg[[m_s_]], title = m_s_, direction = "horizontal")
  }
  if(is.vector(col_leg[[m_s_]])){
    lgd[[m_s_]] <- Legend(labels = names(col_leg[[m_s_]]), legend_gp = gpar(fill = col_leg[[m_s_]]), nrow = 2, title = m_s_)
  }
}

rm(list = ls()[grepl("_$", ls())])
