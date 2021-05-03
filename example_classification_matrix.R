library(corrplot)
cols = colorRampPalette(c('#0000FF', '#FF0000'))
perf = read.csv("example_outcomes/classification_example - Sheet1.csv")
perf = perf[,-1]
x = seq(5, 50, by = 5)
colnames(perf) = x
rownames(perf) = x
perf[is.na(perf)] = 0
corrplot(as.matrix(perf), is.corr = T, method = "square", tl.srt = 0,
         tl.col = "black", cl.pos = "n", tl.offset = 1)


early = read.csv("example_outcomes/classification_example - Sheet2.csv")
early = early[,-1]
colnames(early) = x
rownames(early) = x
early[is.na(early)] = 0


old = read.csv("example_outcomes/classification_example - Sheet3 (1).csv")
old = old[,-1]
colnames(old) = x
rownames(old) = x
old[is.na(old)] = 0

null = read.csv("example_outcomes/classification_example - Sheet4.csv")
null = null[,-1]
colnames(null) = x
rownames(null) = x
null[is.na(null)] = 0

par(mfrow = c(4,1))
corrplot(as.matrix(null), is.corr = T, method = "square", tl.srt = 0,
         tl.col = "black", cl.pos = "n", tl.offset = 1, tl.pos = "n")

corrplot(as.matrix(perf), is.corr = T, method = "square", tl.srt = 0,
         tl.col = "black", cl.pos = "n", tl.offset = 1, tl.pos = "n")

corrplot(as.matrix(early), is.corr = T, method = "square", tl.srt = 0,
         tl.col = "black", cl.pos = "n", tl.offset = 1, tl.pos = "n")

corrplot(as.matrix(old), is.corr = T, method = "square", tl.srt = 0,
         tl.col = "black", cl.pos = "n", tl.offset = 1, tl.pos = "n")

