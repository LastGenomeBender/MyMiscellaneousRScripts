tryCatch({
	library(survival)
	}, error=function(e){
		install.packages("survival",repos="http://cran.pau.edu.tr/")
		library(survival)
	}
)
tryCatch({
	library(maxstat)
	}, error=function(e){
		install.packages("maxstat",repos="http://cran.pau.edu.tr/")
		library(maxstat)
	}
)

setwd("Expression Files")
inputFiles=try(shell("dir /B ",intern=T,wait=T))
setwd("../")

for(files in inputFiles){
	outputUni = paste(substring(files, 1, nchar(files)-4), " Univariate Analysis Results.txt", sep="")
	outputMulti = paste(substring(files, 1, nchar(files)-4), " Multivariate Analysis Results.txt", sep="")
	
	setwd("Expression Files")
	d = as.matrix(read.delim(file=files,sep="\t",header=T))
	setwd("../")
	
	dexpr = d[1:(which(d[,1]=="")[1]-1),]
	dclin = d[(which(d[,1]=="")[1]+1):nrow(d),-ncol(d)]
	dmult1 = c()
	
	for(i in 3:nrow(dclin)){
		dmult1 = cbind(dmult1,as.numeric(dclin[i,-1]))
	}
	colnames(dmult1) = dclin[3:nrow(dclin),1]
	
	surv = as.numeric(dclin[1,-1])
	stat = as.numeric(dclin[2,-1])
	
	resMat = c("Gene","Probeset","HR","Cox p","Threshold","Maxstat p","Median Surv High Expr","Median Surv Low Expr","Log-Rank p","Censored/High Expr","Censored/Low Expr")
	mCoxMat = c("Variables","log(HR)","HR","p-value","Lower 95% Conf","Higher 95% Conf")
	
	for(i in 1:nrow(dexpr)){
		ps = dexpr[i,1]
		gene = dexpr[i,ncol(dexpr)]
		expr = as.numeric(dexpr[i,-c(1,ncol(dexpr))])
		
		hr = ""
		pCox = ""
		
		tryCatch({
			cox = coxph(Surv(surv,stat) ~  expr)
			hr = exp(cox[[1]][[1]])
			pCox = summary(cox)$coefficients[,5][[1]]
		}, error=function(e){})
		
		thr = ""
		pMstat = ""
		pLr = ""
		msLow = ""
		msHigh = ""
		censoredLow = ""
		censoredHigh = ""
		
		tryCatch({
			mstat = maxstat.test(Surv(surv,stat) ~ expr, data=data.frame(expr), smethod="LogRank", pmethod="exactGauss", abseps=0.01)
			thr = mstat$estimate
			pMstat = mstat$p.value
			
			expr2 = expr
			expr2[which(expr<=thr)] = 0
			expr2[which(expr>thr)] = 1
			
			lr = survdiff(Surv(surv,stat) ~  expr2)
			KM = survfit(Surv(surv,stat) ~  expr2)
			msLow = as.vector(summary(KM)$table[,5])[1]
			msHigh = as.vector(summary(KM)$table[,5])[2]
			pLr = pchisq(lr$chisq, 1, lower.tail=FALSE)
			censoredLow = paste(length(which(stat[which(expr<=thr)]==0)),length(which(expr<=thr)), sep="/")
			censoredHigh = paste(length(which(stat[which(expr>thr)]==0)),length(which(expr>thr)), sep="/")
		}, error=function(e){})
		
		resMat = rbind(resMat,c(gene,ps,hr,pCox,thr,pMstat,msHigh,msLow,pLr,censoredHigh,censoredLow))
		
		if(pCox!="" && pCox < 0.05){
			tryCatch({
				dmult2 = cbind(expr,dmult1)
				colnames(dmult2)[1]=paste(gene,ps,sep="_")
				dmult = data.frame(dmult2)
				mCox = coxph(Surv(surv,stat) ~ ., data=dmult)
				mCoxRes = cbind(rownames(summary(mCox)$coef),rbind(summary(mCox)$coef[,c(1,2,5)]),rbind(summary(mCox)$conf[,c(3,4)]))
				mCoxMat = rbind(mCoxMat,"",mCoxRes)
			}, error=function(e){})
		}
	}
	
	setwd("Output Files")
	write(file=outputUni,t(rbind(resMat)),ncol=ncol(rbind(resMat)),sep="\t")
	write(file=outputMulti,t(rbind(mCoxMat)),ncol=ncol(rbind(mCoxMat)),sep="\t")
	setwd("../")
}
