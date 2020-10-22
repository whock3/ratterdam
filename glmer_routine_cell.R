# cellID <- 9
# print(sprintf("Cell ID: %s", cellID))
# mod <- glmer_cell(cellID, df)
# print(sjPlot::plot_model(mod, type='int'))
# print(summary(mod))
# ci = confint(mod, parm='beta_', level = 0.97, method='Wald')
# print("Confidence Intervals (Wald)")
# print(ci)

cellID <- 0
celldf <- subset(df, cell == cellID)
alleys <- unique(celldf$alley)
cellname <- unique(subset(df, cell==cellID, select=name))

for(alleyID in alleys){
  
  print(sprintf("============ Cell: %s, Cell ID: %s, alley ID: %s =============== ", cellname, cellID, alleyID))
  mod <- glmer_alley(cellID, alleyID, df)
  title <- sprintf("Cell %s alley %s", cellname, alleyID)
  print(sjPlot::plot_model(mod, type='int', title=title))
  print(summary(mod))
  ci = confint(mod, parm='beta_', level = 0.97, method='Wald')
  print("Confidence Intervals (Wald)")
  print(ci)
}

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="E:\\Ratterdam\\R781\\glmer_data\\BRD3\\")