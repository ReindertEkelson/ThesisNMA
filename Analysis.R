
DF_RESULTS_ALL <- read.csv("~/Master statistics/DF_RESULTS_ALL.csv")
names(DF_RESULTS_ALL)
attach(DF_RESULTS_ALL)
head(DF_RESULTS_ALL)

CP = DF_RESULTS_ALL[c(1,9,38,29,4,14,19,24,34)]
CPF = cbind(CP[1],round(CP[c(-1)]*100,1))
colnames(CPF) = c("Scenario", "metafor", "brms", "mixmeta", "netmeta", "gemtc U(0,5)", "gemtc gamma(0.001,0.001)", "gemtc HN(1)", "pcnetmeta")

RMSE = DF_RESULTS_ALL[c(1,8,37,28,3,13,18,23,33)]
RMSEF = cbind(RMSE[1],round(RMSE[c(-1)],3))
colnames(RMSEF) = c("Scenario", "metafor", "brms", "mixmeta", "netmeta", "gemtc U(0,5)", "gemtc gamma(0.001,0.001)", "gemtc HN(1)", "pcnetmeta")

MeanBias = DF_RESULTS_ALL[c(1,7,36,27,2,12,17,22,32)]
MeanBiasF = cbind(MeanBias[1],round(MeanBias[c(-1)],4))
colnames(MeanBiasF) = c("Scenario", "metafor", "brms", "mixmeta", "netmeta", "gemtc U(0,5)", "gemtc gamma(0.001,0.001)", "gemtc HN(1)", "pcnetmeta")
###########################

MeanCI = DF_RESULTS_ALL[c(1,10,39,30,5,15,20,25,35)]
MeanCIF = cbind(MeanCI[1],round(MeanCI[c(-1)],3))
colnames(MeanCIF) = c("Scenario", "metafor", "brms", "mixmeta", "netmeta", "gemtc U(0,5)", "gemtc gamma(0.001,0.001)", "gemtc HN(1)", "pcnetmeta")


DATA = expand.grid(data.frame( c(3,5,8), c(3,4,5), c(0.01,(0.1),1)))
colnames(DATA) = c("NS", "NT", "tau")
CPP = cbind(CP, DATA)
CPD = data.frame(a=unlist(CPP[2:9], use.names = FALSE))
Scenario = rbind(CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1])
NS = rbind(CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10])
NT = rbind(CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11])
tau = rbind(CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12])
Estimator = cbind(rep(names(CPP[2:9]), each = 27))
Comb = rep(1:9,24)
PLOTDATA = cbind(CPD, Scenario, NS, NT, tau, Estimator, Comb)
colnames(PLOTDATA) = c("CP", "Scenario", "NS", "NT", "tau", "Estimator", "Comb")


str(PLOTDATA)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
PLOTDATA$Estimator <- factor(PLOTDATA$Estimator, levels = c(names(CPP[2:9])))

CP_Plot = ggplot(PLOTDATA, aes(x = tau, y = CP, colour = Estimator, group = Estimator)) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0.9400, ymax = 0.9600,  alpha = .3) + 
  geom_line(size = 1.5) + geom_point(size = 3) + geom_hline(yintercept = 0.9500, col = "firebrick", lty = 2) + 
  scale_x_log10(breaks=c(0.01,0.1,1)) + 
  scale_color_brewer(palette = "Dark2", name = "Legend",
  label = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U",
            "gemtc gamma", "gemtc HN", "pcnetmeta")) + 
  labs(y = "Coverage Probability of the 95% CI", x = "Between-study variance") + theme_bw() + 
  facet_wrap(~ Comb, labeller = labeller(Comb = c("1" = "NS = 3, NT = 3","2" = "NS = 5, NT = 3","3" = "NS = 8, NT = 3", 
                                                  "4" ="NS = 3, NT = 4","5" = "NS = 5, NT = 4","6" = "NS = 8, NT = 4", 
                                                  "7" ="NS = 3, NT = 5","8" = "NS = 5, NT = 5","9" = "NS = 8, NT = 5"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
                             legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
                             legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) 
CP_Plot


############

RMSE = DF_RESULTS_ALL[c(1,8,37,28,3,13,18,23,33)]
DATA = expand.grid(data.frame( c(3,5,8), c(3,4,5), c(0.01,(0.1),1)))
colnames(DATA) = c("NS", "NT", "tau")
CPP = cbind(RMSE, DATA)
CPD = data.frame(a=unlist(CPP[2:9], use.names = FALSE))
Scenario = rbind(CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1])
NS = rbind(CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10])
NT = rbind(CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11])
tau = rbind(CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12])
Estimator = cbind(rep(names(CPP[2:9]), each = 27))
Comb = rep(1:9,24)
PLOTDATA = cbind(CPD, Scenario, NS, NT, tau, Estimator, Comb)
colnames(PLOTDATA) = c("RMSE", "Scenario", "NS", "NT", "tau", "Estimator", "Comb")
PLOTDATA$Estimator <- factor(PLOTDATA$Estimator, levels = c(names(CPP[2:9])))


RMSE_Plot = ggplot(PLOTDATA, aes(x = tau, y = RMSE, colour = Estimator, group = Estimator)) + 
  geom_line(size = 1.5) + geom_point(size = 3) +
  scale_x_log10(breaks=c(0.01,0.1,1)) + 
  scale_color_brewer(palette = "Dark2", name = "Legend",
                     label = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U",
                               "gemtc gamma", "gemtc HN", "pcnetmeta")) + 
  labs(y = "Root mean squared error (RMSE)", x = "Between-study variance") + theme_bw() + 
  facet_wrap(~ Comb, labeller = labeller(Comb = c("1" = "NS = 3, NT = 3","2" = "NS = 5, NT = 3","3" = "NS = 8, NT = 3", 
                                                  "4" ="NS = 3, NT = 4","5" = "NS = 5, NT = 4","6" = "NS = 8, NT = 4", 
                                                  "7" ="NS = 3, NT = 5","8" = "NS = 5, NT = 5","9" = "NS = 8, NT = 5"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
        legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) 
RMSE_Plot

##################
MeanCI = DF_RESULTS_ALL[c(1,10,39,30,5,15,20,25,35)]
DATA = expand.grid(data.frame( c(3,5,8), c(3,4,5), c(0.01,(0.1),1)))
colnames(DATA) = c("NS", "NT", "tau")
CPP = cbind(MeanCI, DATA)
CPD = data.frame(a=unlist(CPP[2:9], use.names = FALSE))
Scenario = rbind(CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1])
NS = rbind(CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10])
NT = rbind(CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11])
tau = rbind(CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12])
Estimator = cbind(rep(names(CPP[2:9]), each = 27))
Comb = rep(1:9,24)
PLOTDATA = cbind(CPD, Scenario, NS, NT, tau, Estimator, Comb)
colnames(PLOTDATA) = c("MeanCI", "Scenario", "NS", "NT", "tau", "Estimator", "Comb")
PLOTDATA$Estimator <- factor(PLOTDATA$Estimator, levels = c(names(CPP[2:9])))


MeanCI_Plot = ggplot(PLOTDATA, aes(x = tau, y = MeanCI, colour = Estimator, group = Estimator)) + 
  geom_line(size = 1.5) + geom_point(size = 3) +
  scale_x_log10(breaks=c(0.01,0.1,1)) + 
  scale_color_brewer(palette = "Dark2", name = "Legend",
                     label = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U",
                               "gemtc gamma", "gemtc HN", "pcnetmeta")) + 
  labs(y = "Mean width of the 95% CI", x = "Between-study variance") + theme_bw() + 
  facet_wrap(~ Comb, labeller = labeller(Comb = c("1" = "NS = 3, NT = 3","2" = "NS = 5, NT = 3","3" = "NS = 8, NT = 3", 
                                                  "4" ="NS = 3, NT = 4","5" = "NS = 5, NT = 4","6" = "NS = 8, NT = 4", 
                                                  "7" ="NS = 3, NT = 5","8" = "NS = 5, NT = 5","9" = "NS = 8, NT = 5"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
        legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) 
MeanCI_Plot


#################
MeanBias = DF_RESULTS_ALL[c(1,7,36,27,2,12,17,22,32)]
DATA = expand.grid(data.frame( c(3,5,8), c(3,4,5), c(0.01,(0.1),1)))
colnames(DATA) = c("NS", "NT", "tau")
CPP = cbind(MeanBias, DATA)
CPD = data.frame(a=unlist(CPP[2:9], use.names = FALSE))
Scenario = rbind(CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1])
NS = rbind(CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10])
NT = rbind(CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11])
tau = rbind(CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12])
Estimator = cbind(rep(names(CPP[2:9]), each = 27))
Comb = rep(1:9,24)
PLOTDATA = cbind(CPD, Scenario, NS, NT, tau, Estimator, Comb)
colnames(PLOTDATA) = c("MeanBias", "Scenario", "NS", "NT", "tau", "Estimator", "Comb")
PLOTDATA$Estimator <- factor(PLOTDATA$Estimator, levels = c(names(CPP[2:9])))


MeanBias_Plot = ggplot(PLOTDATA, aes(x = tau, y = MeanBias, colour = Estimator, group = Estimator)) + 
  geom_line(size = 1.5) + geom_point(size = 3) +
  scale_x_log10(breaks=c(0.01,0.1,1)) + 
  scale_color_brewer(palette = "Dark2", name = "Legend",
                     label = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U",
                               "gemtc gamma", "gemtc HN", "pcnetmeta")) + 
  labs(y = "Mean Bias", x = "Between-study variance") + theme_bw() + 
  facet_wrap(~ Comb, labeller = labeller(Comb = c("1" = "NS = 3, NT = 3","2" = "NS = 5, NT = 3","3" = "NS = 8, NT = 3", 
                                                  "4" ="NS = 3, NT = 4","5" = "NS = 5, NT = 4","6" = "NS = 8, NT = 4", 
                                                  "7" ="NS = 3, NT = 5","8" = "NS = 5, NT = 5","9" = "NS = 8, NT = 5"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
        legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) 
MeanBias_Plot
########################

tau2 = DF_RESULTS_ALL[c(1,11,40,31,6,16,21,26)] #pcnetmeta could not find a way to extra the tau value
tau2biass = rbind((tau2[1:9,2:8]-0.01), (tau2[10:18,2:8]-0.1), (tau2[19:27,2:8]-1))
tau2bias = cbind(tau2[,1], tau2biass)

tau2biasF = cbind(round(tau2bias[c(-1)],4))
colnames(tau2biasF) = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U(0,5)", "gemtc gamma(0.001,0.001)", "gemtc HN(1)")


DATA = expand.grid(data.frame( c(3,5,8), c(3,4,5), c(0.01,(0.1),1)))
colnames(DATA) = c("NS", "NT", "tau")
CPP = cbind(tau2bias, DATA)
CPD = data.frame(a=unlist(CPP[2:9], use.names = FALSE))
Scenario = rbind(CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1], CPP[1])
NS = rbind(CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10], CPP[10])
NT = rbind(CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11], CPP[11])
tau = rbind(CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12], CPP[12])
Estimator = cbind(rep(names(CPP[2:9]), each = 27))
Comb = rep(1:9,24)
PLOTDATA = cbind(CPD, Scenario, NS, NT, tau, Estimator, Comb)
colnames(PLOTDATA) = c("tau2bias", "Scenario", "NS", "NT", "tau", "Estimator", "Comb")
PLOTDATA$Estimator <- factor(PLOTDATA$Estimator, levels = c(names(CPP[2:9])))
PLOTDATA = PLOTDATA[1:189,]
PLOTDATA[PLOTDATA$tau == 0.01,]

Tau2Bias_Plot1 = ggplot(PLOTDATA[PLOTDATA$tau == 0.01,], aes(x = Estimator, y = tau2bias, fill = Estimator)) + 
  geom_bar(stat='identity', colour = "black") + 
  scale_fill_brewer(palette = "Dark2", name = "Legend",
                     label = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U",
                               "gemtc gamma", "gemtc HN")) + 
  scale_x_discrete(labels = c("metafor", "brms", "mixmeta", "metmeta", "gemtcU", "gemtcG", "gemtcHN")) + 
  labs(y = "Bias of the between-study variance", x = "Between-study variance") + theme_bw() + 
  facet_wrap(~ Comb, labeller = labeller(Comb = c("1" = "NS = 3, NT = 3","2" = "NS = 5, NT = 3","3" = "NS = 8, NT = 3", 
                                                  "4" ="NS = 3, NT = 4","5" = "NS = 5, NT = 4","6" = "NS = 8, NT = 4", 
                                                  "7" ="NS = 3, NT = 5","8" = "NS = 5, NT = 5","9" = "NS = 8, NT = 5"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) 


Tau2Bias_Plot1



Tau2Bias_Plot2 = ggplot(PLOTDATA[PLOTDATA$tau == 0.1,], aes(x = Estimator, y = tau2bias, fill = Estimator)) + 
  geom_bar(stat='identity', colour = "black") + 
  scale_fill_brewer(palette = "Dark2", name = "Legend",
                    label = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U",
                              "gemtc gamma", "gemtc HN")) + 
  scale_x_discrete(labels = c("metafor", "brms", "mixmeta", "netmeta", "gemtcU", "gemtcG", "gemtcHN")) +
  labs(y = "Bias of the between-study variance", x = "Between-study variance") + theme_bw() + 
  facet_wrap(~ Comb, labeller = labeller(Comb = c("1" = "NS = 3, NT = 3","2" = "NS = 5, NT = 3","3" = "NS = 8, NT = 3", 
                                                  "4" ="NS = 3, NT = 4","5" = "NS = 5, NT = 4","6" = "NS = 8, NT = 4", 
                                                  "7" ="NS = 3, NT = 5","8" = "NS = 5, NT = 5","9" = "NS = 8, NT = 5"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) 
Tau2Bias_Plot2




Tau2Bias_Plot3 = ggplot(na.omit(PLOTDATA[PLOTDATA$tau == 1,]), aes(x = Estimator, y = tau2bias, fill = Estimator)) + 
  geom_bar(stat='identity', colour = "black") + 
  scale_fill_brewer(palette = "Dark2", name = "Legend",
                    label = c("netmeta", "gemtc U",
                              "gemtc gamma", "gemtc HN")) + 
  scale_x_discrete(labels = c("metafor", "brms", "mixmeta", "netmeta", "gemtcU", "gemtcG", "gemtcHN")) +
  labs(y = "Bias of the between-study variance", x = "Between-study variance") + theme_bw() + 
  facet_wrap(~ Comb, labeller = labeller(Comb = c("1" = "NS = 3, NT = 3","2" = "NS = 5, NT = 3","3" = "NS = 8, NT = 3", 
                                                  "4" ="NS = 3, NT = 4","5" = "NS = 5, NT = 4","6" = "NS = 8, NT = 4", 
                                                  "7" ="NS = 3, NT = 5","8" = "NS = 5, NT = 5","9" = "NS = 8, NT = 5"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) 
Tau2Bias_Plot3



detach(DF_RESULTS_ALL)

##########################################

DF_Estimates_ALL <- read.csv("~/Master statistics/DF_Estimates_ALL.csv")
attach(DF_Estimates_ALL)
DATA = expand.grid(data.frame( c(3,5,8), c(3,4,5), c(0.01,(0.1),1)))
colnames(DATA) = c("NS", "NT", "tau")
DF_RESULTS_ALLED = cbind(DF_Estimates_ALL, DATA)
DF_RESULTS_ALLED3 = DF_RESULTS_ALLED[DF_RESULTS_ALLED$NT == 3,]
DF_RESULTS_ALLED3 = DF_RESULTS_ALLED3[c(1:3,6:7,10:11,14:15,18:19,22:24,27:28,31:32,35:37)]
DF_RESULTS_ALLED3$X2brms = DF_RESULTS_ALLED3$X1brms + DF_RESULTS_ALLED3$X2brms
DF_RESULTS_ALLED3$X1brms = DF_RESULTS_ALLED3$X1brms + DF_RESULTS_ALLED3$X3brms
DF_RESULTS_ALLED3 = DF_RESULTS_ALLED3[-c(14)]
X = data.frame(matrix(rep(rep(seq(0.5,1,0.5),8), 9), nrow = 9, ncol = 16, byrow = TRUE))
XX = (DF_RESULTS_ALLED3[2:17] - X)

XXX = data.frame(a=unlist(XX))
A = XXX$a
B = rep(DF_RESULTS_ALLED3$tau, 16)
C = rep(DF_RESULTS_ALLED3$NS, 16)
D = cbind(rep(names(XX), each = 9))
E = rep(1:9,16)
F = rep(1:8, each = 18)

FINAL3 = cbind.data.frame(A, B, C, D, E, F)
FINAL3
FINAL3$D <- factor(FINAL3$D, levels = c("X1metafor", "X1brms", "X1mixmeta", "X1netmeta", "X1gemtcU", "X1gemtcG", "X1gemtcHN","X1pcnetmeta",
                                        "X2metafor", "X2brms", "X2mixmeta", "X2netmeta", "X2gemtcU", "X2gemtcG", "X2gemtcHN", "X2pcnetmeta"))

EstimatorPlot3 = ggplot(FINAL3, aes(x = D, y = A, fill = factor(F, levels = c(5,6,7,1,3,2,4,8)))) + 
  geom_bar(stat='identity', colour = "black") +
  scale_fill_brewer(palette = "Dark2", name = "Legend",
                    label = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U",
                              "gemtc gamma", "gemtc HN", "pcnetmeta"))+
  labs(y = "Mean bias", x = "Estimator") + theme_bw() + 
  facet_wrap(~ E, labeller = labeller(E = c("1" = "NS = 3, tau2 = 0.01","2" = "NS = 5, tau2 = 0.01","3" = "NS = 8, tau2 = 0.01", 
                                                  "4" ="NS = 3, tau2 = 0.1","5" = "NS = 5, tau2 = 0.1","6" = "NS = 8, tau2 = 0.1", 
                                                  "7" ="NS = 3, tau2 = 1","8" = "NS = 5, tau2 = 1","9" = "NS = 8, tau2 = 1"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) +
  annotate("text", x = "X1netmeta", y = -0.027, label = 'bold("12")', size = 3, parse = TRUE)+
  annotate("text", x = "X2netmeta", y = -0.027, label = 'bold("13")', size = 3, parse = TRUE)+
  annotate("segment", x = "X1metafor", xend = "X1pcnetmeta", y = -0.025, yend = -0.025, size = 1) + 
  annotate("segment", x = "X2metafor", xend = "X2pcnetmeta", y = -0.025, yend = -0.025, size = 1) + 
  geom_vline(xintercept  = 8.5, col = "firebrick", lty = 2, lwd = 1) 
EstimatorPlot3




DATA = expand.grid(data.frame( c(3,5,8), c(3,4,5), c(0.01,(0.1),1)))
colnames(DATA) = c("NS", "NT", "tau")
DF_RESULTS_ALLED = cbind(DF_Estimates_ALL, DATA)
DF_RESULTS_ALLED4 = DF_RESULTS_ALLED[DF_RESULTS_ALLED$NT == 4,]
DF_RESULTS_ALLED4 = DF_RESULTS_ALLED4[c(1:4,6:8,10:12,14:16,18:20,22:25,27:29,31:33,35:37)]
DF_RESULTS_ALLED4$X3brms = DF_RESULTS_ALLED4$X1brms + DF_RESULTS_ALLED4$X3brms
DF_RESULTS_ALLED4$X2brms = DF_RESULTS_ALLED4$X1brms + DF_RESULTS_ALLED4$X2brms
DF_RESULTS_ALLED4$X1brms = DF_RESULTS_ALLED4$X1brms + DF_RESULTS_ALLED4$X4brms

DF_RESULTS_ALLED4 = DF_RESULTS_ALLED4[-c(20)]
X = data.frame(matrix(rep(rep(seq(1/3,1,1/3),8), 9), nrow = 9, ncol = 24, byrow = TRUE))
XX = (DF_RESULTS_ALLED4[2:25] - X)
XXX = data.frame(a=unlist(XX))
A = XXX$a
B = rep(DF_RESULTS_ALLED4$tau, 24)
C = rep(DF_RESULTS_ALLED4$NS, 24)
D = cbind(rep(names(XX), each = 9))
E = rep(1:9,24)
F = rep(1:8, each = 27)
8*27
FINAL4 = cbind.data.frame(A, B, C, D, E, F)
FINAL4
FINAL4$D <- factor(FINAL4$D, levels = c("X1metafor", "X1brms", "X1mixmeta", "X1netmeta", "X1gemtcU", "X1gemtcG", "X1gemtcHN","X1pcnetmeta",
                                        "X2metafor", "X2brms", "X2mixmeta", "X2netmeta", "X2gemtcU", "X2gemtcG", "X2gemtcHN", "X2pcnetmeta",
                                        "X3metafor", "X3brms", "X3mixmeta", "X3netmeta", "X3gemtcU", "X3gemtcG", "X3gemtcHN", "X3pcnetmeta"))

EstimatorPlot4 = ggplot(FINAL4, aes(x = D, y = A, fill = factor(F, levels = c(5,6,7,1,3,2,4,8)))) + 
  geom_bar(stat='identity', colour = "black") +
  scale_fill_brewer(palette = "Dark2", name = "Legend",
                    label = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U",
                              "gemtc gamma", "gemtc HN", "pcnetmeta"))+
  labs(y = "Mean bias", x = "Estimator") + theme_bw() + 
  facet_wrap(~ E, labeller = labeller(E = c("1" = "NS = 3, tau2 = 0.01","2" = "NS = 5, tau2 = 0.01","3" = "NS = 8, tau2 = 0.01", 
                                            "4" ="NS = 3, tau2 = 0.1","5" = "NS = 5, tau2 = 0.1","6" = "NS = 8, tau2 = 0.1", 
                                            "7" ="NS = 3, tau2 = 1","8" = "NS = 5, tau2 = 1","9" = "NS = 8, tau2 = 1"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) +
  annotate("text", x = "X1netmeta", y = -0.027, label = 'bold("12")', size = 3, parse = TRUE)+
  annotate("text", x = "X2netmeta", y = -0.027, label = 'bold("13")', size = 3, parse = TRUE)+
  annotate("text", x = "X3netmeta", y = -0.027, label = 'bold("14")', size = 3, parse = TRUE)+
  annotate("segment", x = "X1metafor", xend = "X1pcnetmeta", y = -0.025, yend = -0.025, size = 1) + 
  annotate("segment", x = "X2metafor", xend = "X2pcnetmeta", y = -0.025, yend = -0.025, size = 1) + 
  annotate("segment", x = "X3metafor", xend = "X3pcnetmeta", y = -0.025, yend = -0.025, size = 1) + 
  geom_vline(xintercept  = 8.5, col = "firebrick", lty = 2, lwd = 1) +
  geom_vline(xintercept  = 16.5, col = "firebrick", lty = 2, lwd = 1) 
EstimatorPlot4




DATA = expand.grid(data.frame( c(3,5,8), c(3,4,5), c(0.01,(0.1),1)))
colnames(DATA) = c("NS", "NT", "tau")
DF_RESULTS_ALLED = cbind(DF_Estimates_ALL, DATA)
DF_RESULTS_ALLED5 = DF_RESULTS_ALLED[DF_RESULTS_ALLED$NT == 5,]
DF_RESULTS_ALLED5$X4brms = DF_RESULTS_ALLED5$X1brms + DF_RESULTS_ALLED5$X4brms
DF_RESULTS_ALLED5$X3brms = DF_RESULTS_ALLED5$X1brms + DF_RESULTS_ALLED5$X3brms
DF_RESULTS_ALLED5$X2brms = DF_RESULTS_ALLED5$X1brms + DF_RESULTS_ALLED5$X2brms
DF_RESULTS_ALLED5$X1brms = DF_RESULTS_ALLED5$X1brms + DF_RESULTS_ALLED5$X5brms

DF_RESULTS_ALLED5 = DF_RESULTS_ALLED5[-c(26)]
X = data.frame(matrix(rep(rep(seq(0.25,1,0.25),8), 9), nrow = 9, ncol = 32, byrow = TRUE))
XX = (DF_RESULTS_ALLED5[2:33] - X)
XX
XXX = data.frame(a=unlist(XX))
A = XXX$a
B = rep(DF_RESULTS_ALLED5$tau, 32)
C = rep(DF_RESULTS_ALLED5$NS, 32)
D = cbind(rep(names(XX), each = 9))
E = rep(1:9,32)
F = rep(1:8, each = 36)

FINAL5 = cbind.data.frame(A, B, C, D, E, F)
FINAL5
FINAL5$D <- factor(FINAL5$D, levels = c("X1metafor", "X1brms", "X1mixmeta", "X1netmeta", "X1gemtcU", "X1gemtcG", "X1gemtcHN","X1pcnetmeta",
                                        "X2metafor", "X2brms", "X2mixmeta", "X2netmeta", "X2gemtcU", "X2gemtcG", "X2gemtcHN", "X2pcnetmeta",
                                        "X3metafor", "X3brms", "X3mixmeta", "X3netmeta", "X3gemtcU", "X3gemtcG", "X3gemtcHN", "X3pcnetmeta",
                                        "X4metafor", "X4brms", "X4mixmeta", "X4netmeta", "X4gemtcU", "X4gemtcG", "X4gemtcHN", "X4pcnetmeta"))

EstimatorPlot5 = ggplot(FINAL5, aes(x = D, y = A, fill = factor(F, levels = c(5,6,7,1,3,2,4,8)))) + 
  geom_bar(stat='identity', colour = "black") +
  scale_fill_brewer(palette = "Dark2", name = "Legend",
                    label = c("metafor", "brms", "mixmeta", "netmeta", "gemtc U",
                              "gemtc gamma", "gemtc HN", "pcnetmeta"))+
  labs(y = "Mean bias", x = "Estimator") + theme_bw() + 
  facet_wrap(~ E, labeller = labeller(E = c("1" = "NS = 3, tau2 = 0.01","2" = "NS = 5, tau2 = 0.01","3" = "NS = 8, tau2 = 0.01", 
                                            "4" ="NS = 3, tau2 = 0.1","5" = "NS = 5, tau2 = 0.1","6" = "NS = 8, tau2 = 0.1", 
                                            "7" ="NS = 3, tau2 = 1","8" = "NS = 5, tau2 = 1","9" = "NS = 8, tau2 = 1"))) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.width=unit(2,"line"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 11), strip.text.x = element_text(size = 10, face = "bold")) +
  annotate("text", x = "X1netmeta", y = -0.027, label = 'bold("12")', size = 3, parse = TRUE)+
  annotate("text", x = "X2netmeta", y = -0.027, label = 'bold("13")', size = 3, parse = TRUE)+
  annotate("text", x = "X3netmeta", y = -0.027, label = 'bold("14")', size = 3, parse = TRUE)+
  annotate("text", x = "X4netmeta", y = -0.027, label = 'bold("15")', size = 3, parse = TRUE)+
  annotate("segment", x = "X1metafor", xend = "X1pcnetmeta", y = -0.025, yend = -0.025, size = 1) + 
  annotate("segment", x = "X2metafor", xend = "X2pcnetmeta", y = -0.025, yend = -0.025, size = 1) + 
  annotate("segment", x = "X3metafor", xend = "X3pcnetmeta", y = -0.025, yend = -0.025, size = 1) + 
  annotate("segment", x = "X4metafor", xend = "X4pcnetmeta", y = -0.025, yend = -0.025, size = 1) + 
  geom_vline(xintercept  = 8.5, col = "firebrick", lty = 2, lwd = 1) +
  geom_vline(xintercept  = 16.5, col = "firebrick", lty = 2, lwd = 1) +
  geom_vline(xintercept  = 24.5, col = "firebrick", lty = 2, lwd = 1) 
EstimatorPlot5



CP_Plot
RMSE_Plot
MeanCI_Plot
MeanBias_Plot
Tau2Bias_Plot1
Tau2Bias_Plot2
Tau2Bias_Plot3

EstimatorPlot3
EstimatorPlot4
EstimatorPlot5
detach(DF_Estimates_ALL)
