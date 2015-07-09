
# coding: utf-8

# ### Same as Untitled7 but using ranked expr levels instead of actual values

# In[2]:

from matplotlib import pyplot as plt

get_ipython().magic(u'matplotlib inline')
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('pdf', 'svg')

import pandas as pd
pd.options.display.mpl_style = 'default'

from mpltools import style
from mpltools import layout

style.use('ggplot')

## see: https://stackoverflow.com/questions/19536817/manipulate-html-module-font-size-in-ipython-notebook
class sizeme():
    """ Class to change html fontsize of object's representation"""
    def __init__(self,ob, size, height=100):
        self.ob = ob
        self.size = size
        self.height = height
    def _repr_html_(self):
        repl_tuple = (self.size, self.height, self.ob._repr_html_())
        return u'<span style="font-size:{0}%; line-height:{1}%">{2}</span>'.format(*repl_tuple)
    
## see https://stackoverflow.com/questions/14656852/how-to-use-pandas-dataframes-and-numpy-arrays-in-rpy2
## and http://ipython.org/ipython-doc/rel-0.13/config/extensions/rmagic.html
## note there's a ri2pandas() to convert back.
## but note, rpy2 2.4.0 and later automagically translates dataframes: 
## https://stackoverflow.com/questions/20630121/pandas-how-to-convert-r-dataframe-back-to-pandas
get_ipython().magic(u'load_ext rpy2.ipython')
get_ipython().magic(u'Rdevice svg')
#import rpy2.robjects.pandas2ri as p2r
#rdf = p2r.pandas2ri(info)
#%Rpush rdf


# In[ ]:

get_ipython().magic(u'run -i read_counts.py')


# In[31]:

sample_info = pd.read_excel('Sample_Info_FIXED2.xlsx') ##,skiprows=[0])
sample_info = sample_info.drop( ['growth rate per h (OLD)', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Sample name.1'], 1)
sample_info = sample_info.set_index( sample_info['Sample name'] )
sample_infos = { k:sample_info.ix[all_freqs[k].columns.droplevel(1).values] for k in all_freqs.keys() }
info = sample_infos['Desulfovibrio_vulgaris_Hildenborough_uid57645'].copy()
sizeme(info.head(3),50,120)


# ### Identify groups of replicates in the measurements -- using groupby

# In[33]:

#info_tmp = info[info.columns[np.hstack([4,np.arange(6,14)])]]  ##.duplicated()
group_cols = info.columns[np.hstack([4,6,7,8,9,11,12,13,14])].values.astype(str).tolist()
#grouped = info.groupby(info.columns[np.hstack([4,np.arange(6,14)])].values)
#group_cols = ['cultivation type', 'Description/condition details',
#             'Description/condition details -2', 'carbon source', 'electron donor',
#             'concentration (mM)', 'electron acceptor', 'growth rate per h', 'T0C', 'organisms']
print group_cols
grouped = info.groupby(group_cols, axis=0)
print len(grouped.groups), info.shape, info[group_cols].drop_duplicates().shape
#print grouped.groups[grouped.groups.keys()[0]]
print [len(i) for i in grouped.groups.values()]
col_groups = grouped.groups.values()
for i in grouped.groups.values():
    print i
    #print info.ix[i]['Description/condition details'].values
col_groups = grouped.groups.values()


# ## OK, idea: remove all replicates of a given measurement, run Boruta to get big, best subset of variables that classify, then random forest (lots of trees) using those variables to get classifier. Test (predict) the random forest on the left-out measurements.

# ### Now do it for all test cases! -- i.e., all replicate sets

# In[6]:

get_ipython().run_cell_magic(u'R', u'-n', u"load('Untitled5.RData')\nload('Untitled6.RData')\ncolnames(x) <- gsub('X..','',gsub('...1.','',colnames(x),fixed=T),fixed=T)\nrequire(Boruta); require(randomForest); require(parallel); options(mc.cores=8); options(cores=8)\n\ndo_it <- function(cond_type, x, cols_exclude=NULL, genes_exclude=NULL, n_trees=100000) {\n    if (cond_type == 'growth_rate') {\n        gr <- info$growth.rate.per.h; gr=as.numeric(as.character(gr)); gr[is.na(gr)]=0.2\n        Y <- factor( ifelse(gr==0, 'no_growth', ifelse(gr>=0.2, 'fast_growth', 'med_growth')) )\n    } else if (cond_type == 'electron_donor') {\n        Y <- as.factor(as.character(info$electron.donor))\n    } else if (cond_type == 'electron_acceptor') {\n        Y <- as.factor(as.character(info$electron.acceptor))\n    } else if (cond_type == 'temperature') {\n        Y <- as.factor(as.character(info$T0C == 37))\n    } else if (cond_type == 'all') {\n        Y1 <- '' #ifelse(gr==0, 'no_growth', ifelse(gr>=0.2, 'fast_growth', 'med_growth'))\n        Y2 <- as.character(info$electron.donor)\n        Y3 <- as.character(info$electron.acceptor)\n        Y4 <- as.character(info$T0C == 37)\n        Y <- as.factor(paste(Y1, Y2, Y3, Y4))\n    }\n    #cat(cond_type, length(levels(Y)), '\\n')\n    names(Y) <- info$Sample.name\n\n    cols2 <- ''\n    if ( ! is.null(cols_exclude) ) cols2 <- gsub('-','.',cols_exclude,fixed=T)\n\n    YY <- Y\n    if ( ! is.null(cols_exclude) ) YY <- Y[!names(Y)%in%cols_exclude]\n        \n    XX <- x\n    if ( ! is.null(cols_exclude) ) XX <- XX[, !colnames(XX) %in% cols2, drop=F]\n    if ( ! is.null(genes_exclude) ) XX <- XX[!rownames(XX) %in% genes_exclude, ,drop=F]\n    #print(dim(XX));print(length(YY))\n\n    B.temp1a <- Boruta(t(XX), YY, getImp=getImpFerns, ferns=n_trees, doTrace=0)\n    features <- gsub('.','-',getSelectedAttributes(B.temp1a), fixed=T)\n        \n    rf.temp1a <- randomForest(t(XX[features, ,drop=F]), YY, importance=T, ntree=n_trees, do.trace=F)\n    tmp <- list(predicted=predict(rf.temp1a, t(x[features,cols2,drop=F])), actual=Y[cols_exclude], features=features)\n    \n    return(tmp)\n}\nNULL")


# In[35]:

get_ipython().magic(u'R -n x.ranks <- apply(x,2,rank)')
get_ipython().magic(u'Rpush col_groups')
get_ipython().magic(u'R print(length(col_groups))')
#%R save(col_groups,file='qqq')


# In[9]:

get_ipython().run_cell_magic(u'R', u'', u"results = list()\nfor (cond_type in c('electron_donor', 'electron_acceptor', 'temperature', 'growth_rate')) {\n    tmp <- mclapply( col_groups, function(cols) {\n        cols = unlist(cols)\n        tmp <- do_it(cond_type, x.ranks, cols)\n        #cat(cond_type, mean(as.character(tmp$predicted) == as.character(tmp$actual)), '\\n')\n        return(tmp)\n    }, mc.preschedule=F )\n    cat(cond_type, mean(unlist(lapply(tmp,'[[','predicted')) == unlist(lapply(tmp,'[[','actual'))), '\\n')\n    results[[cond_type]] <- tmp\n}")


# In[10]:

get_ipython().magic(u"R print(sapply(results,function(tmp)mean(unlist(lapply(tmp,'[[','predicted')) == unlist(lapply(tmp,'[[','actual')))))")
get_ipython().magic(u"R print(sapply(results,function(tmp)length(levels(unlist(lapply(tmp,'[[','actual'))))))")
## get table of how many times a given feature was chosen out of the 24 condition replicate "groups":
## %R lapply(lapply(results,lapply,'[[','features'),function(i)sort(table(unlist(i))))
get_ipython().magic(u'Rpull results')


# ### See which classes are most misclassified (growth_rate).
# #### looks like fast growth is most often misclassified as no growth (and vice versa).

# In[11]:

get_ipython().run_cell_magic(u'R', u'-n', u"tmp1=unlist(lapply(lapply(results$growth_rate,'[[','predicted'),as.character))\ntmp2=unlist(lapply(lapply(results$growth_rate,'[[','actual'),as.character))\nprint(table(tmp2[tmp1!=tmp2]))   ## which actual classes are most frequently misclassified\nprint(table(tmp1[tmp1!=tmp2]))   ## which predicted classes are they most frequently misclassified as\ncbind(predicted=tmp1[tmp1!=tmp2],actual=tmp2[tmp1!=tmp2])   ## misclassification pairs (predicted, actual)")


# ### For comparison, do it for the full data sets in which we did not leave any conditions out

# In[ ]:

get_ipython().run_cell_magic(u'R', u'', u"    results_noleaveout = list()\n    for (cond_type in c('growth_rate', 'electron_donor', 'electron_acceptor', 'temperature')) {\n        if (cond_type == 'growth_rate') {\n            gr <- info$growth.rate.per.h; gr=as.numeric(as.character(gr)); gr[is.na(gr)]=0.2\n            Y <- factor( ifelse(gr==0, 'no_growth', ifelse(gr>=0.2, 'fast_growth', 'med_growth')) )\n        } else if (cond_type == 'electron_donor') {\n            Y <- as.factor(as.character(info$electron.donor))\n        } else if (cond_type == 'electron_acceptor') {\n            Y <- as.factor(as.character(info$electron.acceptor))\n        } else if (cond_type == 'temperature') {\n            Y <- as.factor(as.character(info$T0C == 37))\n        } else if (cond_type == 'all') {\n            Y1 <- '' #ifelse(gr==0, 'no_growth', ifelse(gr>=0.2, 'fast_growth', 'med_growth'))\n            Y2 <- as.character(info$electron.donor)\n            Y3 <- as.character(info$electron.acceptor)\n            Y4 <- as.character(info$T0C == 37)\n            Y <- as.factor(paste(Y1, Y2, Y3, Y4))\n        }\n        cat(cond_type, length(levels(Y)), '\\n')\n\n        B.temp1a <- Boruta(t(x.ranks), Y, getImp=getImpFerns, ferns=100000, doTrace=0)\n        features <- gsub('.','-',getSelectedAttributes(B.temp1a), fixed=T)\n        \n        rf.temp1a <- randomForest(t(x.ranks[features,]), Y, importance=T, ntree=100000, do.trace=F)\n        tmp = list(predicted=predict(rf.temp1a), actual=Y, features=features)\n\n        ##print(lapply(tmp,function(i)mean(i$predicted==i$actual)))\n        cat(cond_type, mean(as.character(tmp$predicted) == as.character(tmp$actual)), '\\n')\n        results_noleaveout[[cond_type]] <- tmp\n    }")


# In[ ]:

get_ipython().magic(u"R print(apply(sapply(results_noleaveout,'[[','predicted') == sapply(results_noleaveout,'[[','actual'),2,mean))")
get_ipython().magic(u'R print(sapply(results_noleaveout,function(tmp)length(levels(tmp$actual))))')
## get table of how many times a given feature was chosen out of the 24 condition replicate "groups":
## %R lapply(lapply(results,lapply,'[[','features'),function(i)sort(table(unlist(i))))
get_ipython().magic(u'Rpull results_noleaveout')
get_ipython().magic(u'R rm(rf.temp1a,B.temp1a); save.image("Untitled7.RData")')


# ### Let's try the removing of the lowest-expressed genes incrementally and see how training/testing works for electron_donor conditions:
# 
# ##### lets do it by quantiles first, just to see (try 8 quantiles???):

# In[ ]:

get_ipython().run_cell_magic(u'R', u'', u"require(Boruta); require(randomForest); require(parallel); options(mc.cores=8); options(cores=8)\n#cond_type <- 'electron_donor'\nqqq <- list()\ninp_quantiles <- c(0,0.5,0.75, seq(0.8,0.95,by=0.5),seq(0.92,0.98,by=0.02),seq(0.982,0.999,by=0.002))\nfor (cond_type in c('electron_donor', 'electron_acceptor', 'temperature', 'growth_rate')) {\ntmp_lst <- lapply( 1:length(col_groups), function(col) {\n    #cols = unlist(col_groups[[13]])\n    cols <- col_groups[[col]]\n    cols <- unlist(cols)\n    cols2 <- gsub('-','.',cols)   ## try a pyruvate condition  -- all predicted correctly\n    mns <- apply(x.ranks[,cols2, drop=F], 1, median)\n    levels <- quantile(mns, inp_quantiles)\n    tmp <- mclapply( rev(levels), function(lev) {\n        cat(col, 'of', length(col_groups), rev(names(levels))[which(rev(levels)==lev)], lev, sum(mns>lev), cond_type, '\\n')\n        xx <- x.ranks[mns > lev,, drop=F]\n        ttmp <- do_it(cond_type, xx, cols)\n        #cat(lev, sum(mns>lev), cond_type, mean(as.character(ttmp$predicted) == as.character(ttmp$actual)), '\\n')\n        return(ttmp)\n    }, mc.preschedule=F )\n    #print(apply(sapply(tmp,'[[','predicted') == sapply(tmp,'[[','actual'),2,mean))\n    return( tmp )\n} )\nprint(sapply(tmp_lst,function(tmp)mean(unlist(lapply(tmp,'[[','predicted')) == unlist(lapply(tmp,'[[','actual')))))\nqqq[[cond_type]] <- tmp_lst\n}")


# In[ ]:

get_ipython().magic(u"R save.image('Untitled8.RData')")
#%R load('Untitled8.RData')


# In[25]:

get_ipython().run_cell_magic(u'R', u'', u"load('Untitled8.RData')\ninp_quantiles <- c(0,0.5,0.75, seq(0.8,0.95,by=0.5),seq(0.92,0.98,by=0.02),seq(0.982,0.999,by=0.002))\nmns <- sapply(1:length(inp_quantiles),function(lev)sapply(qqq,function(tmp1)mean(unlist(sapply(tmp1,function(tmp2)tmp2[[lev]]$predicted==unlist(tmp2[[lev]]$actual))))))\ncolnames(mns)<-names(qqq$growth_rate[[1]])\nmns<-mns[,ncol(mns):1]\nexpr_levels = quantile(apply(x.ranks,1,median),inp_quantiles)\nn_genes <- sapply(expr_levels,function(lev)sum(apply(x.ranks,1,median)>lev))\n\nn_features <- sapply(1:length(inp_quantiles),function(lev)sapply(qqq,function(tmp1)mean(unlist(sapply(tmp1,function(tmp2)length(tmp2[[lev]]$features))))))\ncolnames(n_features)<-names(qqq$growth_rate[[1]])\nn_features<-n_features[,ncol(n_features):1]\n\nmatplot(log10(1-inp_quantiles),t(mns),typ='l',main='accuracy')\nlegend('topleft',legend=rownames(mns),lwd=1,col=1:4)\n#text(log10(1-inp_quantiles), 0.5, lab=as.character(round(expr_levels)),cex=0.5,srt=90, xpd=NA)\n#text(0.1, 0.5, lab='Rank', pos=4, cex=0.5, xpd=NA)\ntext(log10(1-inp_quantiles), 0.53, lab=as.character(n_genes),cex=0.5,srt=90, xpd=NA)\ntext(0.1, 0.53, lab='Ngenes', pos=4, cex=0.5, xpd=NA)\ntext(log10(1-inp_quantiles), 0.56, lab=as.character(round(apply(n_features,2,mean))),cex=0.5,srt=90, xpd=NA)\ntext(0.1, 0.56, lab='Nfeatures', pos=4, cex=0.5, xpd=NA)\n\nNULL")


# In[22]:

get_ipython().run_cell_magic(u'R', u'', u"hist(log10(apply(x.ranks,1,median)+1),breaks=200,main='median ranks')\ntmp<-sapply(expr_levels,function(i)lines(rep(log10(i+1),2),c(-1,0),col='red',lwd=3));rm(tmp);NULL")


# In[ ]:



