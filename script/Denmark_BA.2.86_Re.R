#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(ggplot2)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)
library(scales)


##########args##########
#input

stan_f.name <- 'script/multinomial_independent.stan'
metadata.name <- "metadata_tsv_2023_09_04/metadata.tsv"

##########parameters##########
#general
core.num <- 4

#min numbers
limit.count.analyzed <- 50

#Transmissibility
bin.size <- 1
generation_time <- 2.1
pango.reference <- "XBB.1.5"
date.start <- as.Date("2023-04-01")

#model
multi_nomial_model <- cmdstan_model(stan_f.name)


metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)

metadata.filtered <- metadata %>%
                       distinct(Accession.ID,.keep_all=T) %>%
                       filter(Host == "Human",
                              Collection.date > date.start, 
                              str_length(Collection.date) == 10,
                              Pango.lineage != "",
                              Pango.lineage != "None",
                              Pango.lineage != "Unassigned",
                              !str_detect(Additional.location.information,"[Qq]uarantine")
                              )
 
metadata.filtered <- metadata.filtered %>%
                       mutate(Collection.date = as.Date(Collection.date),
                              region = str_split(Location," / ",simplify = T)[,1],
                              country = str_split(Location," / ",simplify = T)[,2],
                              state = str_split(Location," / ",simplify = T)[,3])

metadata.filtered.interest <- metadata.filtered %>% filter(country == "Denmark")
count_country_pango.df <- metadata.filtered.interest %>% group_by(country,Pango.lineage) %>% summarize(count_country_pango = n()) %>% filter(count_country_pango > limit.count.analyzed)

lineage.analyzed.v <- c(count_country_pango.df$Pango.lineage, "BA.2.86")
lineage.interest.v <- c("XBB.1.5","XBB.1.16","EG.5.1","BA.2.86")



metadata.filtered.interest <- metadata.filtered.interest %>% filter(Pango.lineage %in% as.character(lineage.analyzed.v))
metadata.filtered.interest <- metadata.filtered.interest %>% mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1, date.bin = cut(date.num,seq(0,max(date.num),bin.size)), date.bin.num = as.numeric(date.bin))
metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))
#metadata.filtered.final <- rbind(metadata.filtered.final,metadata.filtered.interest)

metadata.filtered.interest.bin <- metadata.filtered.interest %>% group_by(date.bin.num,Pango.lineage) %>% summarize(count = n()) %>% ungroup()

metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% spread(key=Pango.lineage,value = count)

metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0
metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin.spread

X <- as.matrix(data.frame(X0 = 1, X1 = metadata.filtered.interest.bin.spread$date.bin.num))

Y <- metadata.filtered.interest.bin.spread %>% select(- date.bin.num)
Y <- Y[,c(pango.reference,colnames(Y)[-which(colnames(Y)==pango.reference)])]


count.group <- apply(Y,2,sum)
count.total <- sum(count.group)
prop.group <- count.group / count.total

Y <- Y %>% as.matrix()
apply(Y,2,sum)



group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))

Y_sum.v <- apply(Y,1,sum)


data.stan <- list(K = ncol(Y),
                  D = 2,
                  N = nrow(Y),
                  X = X,
                  Y = Y,
                  generation_time = generation_time,
                  bin_size = bin.size,
                  Y_sum = c(Y_sum.v))

fit.stan <- multi_nomial_model$sample(
  data=data.stan,
  iter_sampling=4000,
  iter_warmup=1000,
  seed=1234,
  parallel_chains = 4,
  max_treedepth = 15,
  chains=4)


#growth rate
stat.info <- fit.stan$summary("growth_rate") %>% as.data.frame()
stat.info$Nextclade_pango <- colnames(Y)[2:ncol(Y)]

stat.info.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame() %>% rename(q2.5 = `2.5%`, q97.5 = `97.5%`)
stat.info <- stat.info %>% inner_join(stat.info.q,by="variable")


draw.df.growth_rate <- fit.stan$draws("growth_rate", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.growth_rate.long <- draw.df.growth_rate %>% gather(key = class, value = value)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group_Id = str_match(draw.df.growth_rate.long$class,'growth_rate\\[([0-9]+)\\]')[,2] %>% as.numeric() + 1)
draw.df.growth_rate.long <- merge(draw.df.growth_rate.long,group.df,by="group_Id") %>% select(value,group)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% group_by(group) %>% filter(value>=quantile(value,0.005),value<=quantile(value,0.995))
draw.df.growth_rate.long <- rbind(data.frame(group=pango.reference,value=1),draw.df.growth_rate.long)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% filter(group %in% lineage.interest.v)

draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group = factor(group,levels=lineage.interest.v))

col.v <- brewer.pal(length(lineage.interest.v), "Set1")
draw.df.growth_rate.long <- draw.df.growth_rate.long
g <- ggplot(draw.df.growth_rate.long,aes(x=group,y=value,color=group,fill=group))
g <- g + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
g <- g + geom_violin(alpha=0.6,scale="width")
g <- g + stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975), size=0.5,fatten =1.5)
g <- g + scale_color_manual(values=col.v)
g <- g + scale_fill_manual(values=col.v)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8))
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + xlab('') + ylab(paste('Relative Re (',pango.reference,')',sep=""))
g <- g + theme(legend.position = 'none')
g

pdf.name <- "BA.2.86_Re.pdf"
pdf(pdf.name,width=2,height=3)
plot(g)
dev.off()


#theta
data_Id.df <- data.frame(data_Id = 1:length(X), date_Id = X[,2], Y_sum = Y_sum.v, date = as.Date(X[,2],origin=date.start)-1)


data.freq <- metadata.filtered.interest.bin %>% rename(group = Pango.lineage) %>% group_by(date.bin.num) %>% mutate(freq = count / sum(count))
data.freq <- data.freq %>% mutate(date = as.Date(date.bin.num,origin=date.start)-1)


draw.df.theta <- fit.stan$draws("theta", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.theta.long <- draw.df.theta %>% gather(key = class, value = value)
draw.df.theta.long <- draw.df.theta.long %>% mutate(data_Id = str_match(class,'theta\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(),
                                                    group_Id = str_match(class,'theta\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric())

draw.df.theta.long <- draw.df.theta.long %>% inner_join(data_Id.df %>% select(data_Id,date), by = "data_Id")



draw.df.theta.long.sum <- draw.df.theta.long %>% group_by(group_Id, date) %>% summarize(mean = mean(value),ymin = quantile(value,0.025),ymax = quantile(value,0.975))
draw.df.theta.long.sum <- draw.df.theta.long.sum %>% inner_join(group.df,by="group_Id")

draw.df.theta.long.sum.filtered <- draw.df.theta.long.sum %>% filter(group %in% lineage.interest.v) %>% mutate(group = factor(group,levels=lineage.interest.v))

g <- ggplot(draw.df.theta.long.sum.filtered,aes(x=date, y = mean, fill=group, color = group))
g <- g + geom_ribbon(aes(ymin=ymin,ymax=ymax), color=NA,alpha=0.4)
g <- g + geom_line(size=0.3)
g <- g + scale_x_date(date_labels = "%y-%m", date_breaks = "1 months", date_minor_breaks = "1 month", limits=c(date.start,as.Date("2023-08-31")))
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + scale_color_manual(values = col.v)
g <- g + scale_fill_manual(values = col.v)

pdf.name <- "BA.2.86_theta.pdf"
pdf(pdf.name,width=4.5,height=3)
plot(g)
dev.off()


out.name <- "Re_BA.2.86.txt"
write.table(stat.info,out.name,col.names=T,row.names=F,sep="\t",quote=F)

out.name <- "metadata_BA.2.86.txt"
write.table(metadata.filtered.interest,out.name,col.names=T,row.names=F,sep="\t",quote=F)


draw.df.growth_rate <- fit.stan$draws("growth_rate", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.growth_rate.long <- draw.df.growth_rate %>% mutate(sampling_Id = 1:nrow(draw.df.growth_rate)) %>% gather(key = class, value = value, - sampling_Id)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group_Id = str_match(draw.df.growth_rate.long$class,'growth_rate\\[([0-9]+)\\]')[,2] %>% as.numeric() + 1)
draw.df.growth_rate.long <- merge(draw.df.growth_rate.long,group.df,by="group_Id") %>% select(value,group,sampling_Id)
draw.df.growth_rate.long.spread <- draw.df.growth_rate.long %>% filter(group %in% c("BA.2.86","EG.5.1")) %>% spread(key=group,value = value)
draw.df.growth_rate.long.spread <- draw.df.growth_rate.long.spread %>% mutate(vs_EG.5.1 = BA.2.86 / EG.5.1)

prop.gt1 <- nrow(draw.df.growth_rate.long.spread %>% filter(vs_EG.5.1 > 1)) / nrow(draw.df.growth_rate.long.spread)
prop.gt1 <- round(prop.gt1,3)

g <-  ggplot(draw.df.growth_rate.long.spread, aes(x=vs_EG.5.1, y = ..density..))
g <- g + geom_histogram(binwidth=0.01, fill="gray70")
g <- g + geom_vline(xintercept=1, color = "brown")
g <- g + annotate("text",x=-Inf,y=Inf,label=as.character(prop.gt1),hjust=-.2,vjust=2)
g

pdf.name <- "BA.2.86_vs_EG.5.1.pdf"
pdf(pdf.name,width=3,height=3)
plot(g)
dev.off()


