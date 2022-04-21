
library(networkD3)
data(MisLinks)
data(MisNodes)

MisLinks$value=rep(3,nrow(MisLinks))






MisNodes$size=rep(10,77)
del_label<-c(which(MisNodes$group==0),which(MisNodes$group==7),which(MisNodes$group==10),which(MisNodes$group==6),
             which(MisNodes$group==9))

del_labelLinks_list<-NULL
for (i in 1:length(del_label)) {
  del_labelLinks<-del_label-1
  del_labelLinks_list<-append(del_labelLinks_list,
                            c(which(MisLinks$source==del_labelLinks[i]),which(MisLinks$target==del_labelLinks[i])))
  
}

del_labelLinks_list=unique(del_labelLinks_list)


MisLinks=MisLinks[-del_labelLinks_list,]


del_2_3<-unique(c(which(MisLinks$source==2), which(MisLinks$source==3), which(MisLinks$target==2), which(MisLinks$target==3)))

MisLinks=MisLinks[-del_2_3,]



MisLinks=rbind(c(1,4,3), MisLinks)
MisLinks=rbind(c(1,5,3), MisLinks)
MisLinks=rbind(c(1,6,3), MisLinks)
MisLinks=rbind(c(1,7,3), MisLinks)
MisLinks=rbind(c(1,8,3), MisLinks)
MisLinks=rbind(c(1,9,3), MisLinks)

MisLinks=rbind(c(4,5,3), MisLinks)
MisLinks=rbind(c(4,6,3), MisLinks)
MisLinks=rbind(c(4,7,3), MisLinks)
MisLinks=rbind(c(4,8,3), MisLinks)
MisLinks=rbind(c(4,9,3), MisLinks)

MisLinks=rbind(c(5,6,3), MisLinks)
MisLinks=rbind(c(5,7,3), MisLinks)
MisLinks=rbind(c(5,8,3), MisLinks)
MisLinks=rbind(c(5,9,3), MisLinks)

MisLinks=rbind(c(6,7,3), MisLinks)
MisLinks=rbind(c(6,8,3), MisLinks)
MisLinks=rbind(c(6,9,3), MisLinks)

MisLinks=rbind(c(7,8,3), MisLinks)
MisLinks=rbind(c(7,9,3), MisLinks)

MisLinks=rbind(c(8,9,3), MisLinks)

MisLinks=rbind(c(1,11,3), MisLinks)
MisLinks=rbind(c(1,10,3), MisLinks)
MisLinks=rbind(c(1,13,3), MisLinks)
MisLinks=rbind(c(0,13,3), MisLinks)


MisLinks=rbind(c(10,11,3), MisLinks)
MisLinks=rbind(c(10,13,3), MisLinks)
MisLinks=rbind(c(10,14,3), MisLinks)
MisLinks=rbind(c(10,15,3), MisLinks)
MisLinks=rbind(c(10,29,3), MisLinks)

MisLinks=rbind(c(11,13,3), MisLinks)
MisLinks=rbind(c(11,14,3), MisLinks)
MisLinks=rbind(c(11,15,3), MisLinks)
MisLinks=rbind(c(11,29,3), MisLinks)

MisLinks=rbind(c(13,14,3), MisLinks)
MisLinks=rbind(c(13,15,3), MisLinks)
MisLinks=rbind(c(13,29,3), MisLinks)

MisLinks=rbind(c(14,15,3), MisLinks)
MisLinks=rbind(c(14,29,3), MisLinks)


MisLinks=rbind(c(10,31,3), MisLinks)
MisLinks=rbind(c(10,32,3), MisLinks)
MisLinks=rbind(c(10,33,3), MisLinks)
MisLinks=rbind(c(10,34,3), MisLinks)
MisLinks=rbind(c(10,35,3), MisLinks)
MisLinks=rbind(c(10,36,3), MisLinks)
MisLinks=rbind(c(10,37,3), MisLinks)
MisLinks=rbind(c(10,38,3), MisLinks)


MisLinks=rbind(c(11,31,3), MisLinks)
MisLinks=rbind(c(11,32,3), MisLinks)
MisLinks=rbind(c(11,33,3), MisLinks)
MisLinks=rbind(c(11,34,3), MisLinks)
MisLinks=rbind(c(11,35,3), MisLinks)
MisLinks=rbind(c(11,36,3), MisLinks)
MisLinks=rbind(c(11,37,3), MisLinks)
MisLinks=rbind(c(11,38,3), MisLinks)

MisLinks=rbind(c(13,31,3), MisLinks)
MisLinks=rbind(c(13,32,3), MisLinks)
MisLinks=rbind(c(13,33,3), MisLinks)
MisLinks=rbind(c(13,34,3), MisLinks)
MisLinks=rbind(c(13,35,3), MisLinks)
MisLinks=rbind(c(13,36,3), MisLinks)
MisLinks=rbind(c(13,37,3), MisLinks)
MisLinks=rbind(c(13,38,3), MisLinks)




MisLinks=rbind(c(12,16,3), MisLinks)
MisLinks=rbind(c(12,17,3), MisLinks)
MisLinks=rbind(c(12,18,3), MisLinks)
MisLinks=rbind(c(12,19,3), MisLinks)
MisLinks=rbind(c(12,20,3), MisLinks)
MisLinks=rbind(c(12,21,3), MisLinks)
MisLinks=rbind(c(12,22,3), MisLinks)
MisLinks=rbind(c(12,23,3), MisLinks)



MisLinks=rbind(c(30,16,3), MisLinks)
MisLinks=rbind(c(30,17,3), MisLinks)
MisLinks=rbind(c(30,18,3), MisLinks)
MisLinks=rbind(c(30,19,3), MisLinks)
MisLinks=rbind(c(30,20,3), MisLinks)
MisLinks=rbind(c(30,21,3), MisLinks)
MisLinks=rbind(c(30,22,3), MisLinks)
MisLinks=rbind(c(30,23,3), MisLinks)


MisLinks=rbind(c(24,41,3), MisLinks)
MisLinks=rbind(c(24,42,3), MisLinks)
MisLinks=rbind(c(24,68,3), MisLinks)
MisLinks=rbind(c(24,69,3), MisLinks)
MisLinks=rbind(c(24,70,3), MisLinks)
MisLinks=rbind(c(24,71,3), MisLinks)
MisLinks=rbind(c(24,76,3), MisLinks)


MisLinks=rbind(c(25,41,3), MisLinks)
MisLinks=rbind(c(25,42,3), MisLinks)
MisLinks=rbind(c(25,68,3), MisLinks)
MisLinks=rbind(c(25,69,3), MisLinks)
MisLinks=rbind(c(25,70,3), MisLinks)
MisLinks=rbind(c(25,71,3), MisLinks)
MisLinks=rbind(c(25,76,3), MisLinks)



MisLinks=rbind(c(27,41,3), MisLinks)
MisLinks=rbind(c(27,42,3), MisLinks)
MisLinks=rbind(c(27,68,3), MisLinks)
MisLinks=rbind(c(27,69,3), MisLinks)
MisLinks=rbind(c(27,70,3), MisLinks)
MisLinks=rbind(c(27,71,3), MisLinks)
MisLinks=rbind(c(27,76,3), MisLinks)

MisLinks=rbind(c(29,41,3), MisLinks)
MisLinks=rbind(c(29,42,3), MisLinks)
MisLinks=rbind(c(29,68,3), MisLinks)
MisLinks=rbind(c(29,69,3), MisLinks)
MisLinks=rbind(c(29,70,3), MisLinks)
MisLinks=rbind(c(29,71,3), MisLinks)
MisLinks=rbind(c(29,76,3), MisLinks)



MisLinks=rbind(c(41,42,3), MisLinks)
MisLinks=rbind(c(41,68,3), MisLinks)
MisLinks=rbind(c(41,69,3), MisLinks)
MisLinks=rbind(c(41,70,3), MisLinks)
MisLinks=rbind(c(41,71,3), MisLinks)
MisLinks=rbind(c(41,76,3), MisLinks)

MisLinks=rbind(c(42,68,3), MisLinks)
MisLinks=rbind(c(42,69,3), MisLinks)
MisLinks=rbind(c(42,70,3), MisLinks)
MisLinks=rbind(c(42,71,3), MisLinks)
MisLinks=rbind(c(42,76,3), MisLinks)

MisLinks=rbind(c(26,49,3), MisLinks)
MisLinks=rbind(c(26,50,3), MisLinks)
MisLinks=rbind(c(26,51,3), MisLinks)
MisLinks=rbind(c(26,52,3), MisLinks)
MisLinks=rbind(c(26,53,3), MisLinks)
MisLinks=rbind(c(26,54,3), MisLinks)


MisLinks=rbind(c(43,49,3), MisLinks)
MisLinks=rbind(c(43,50,3), MisLinks)
MisLinks=rbind(c(43,51,3), MisLinks)
MisLinks=rbind(c(43,52,3), MisLinks)
MisLinks=rbind(c(43,53,3), MisLinks)
MisLinks=rbind(c(43,54,3), MisLinks)



MisLinks=rbind(c(56,49,3), MisLinks)
MisLinks=rbind(c(56,50,3), MisLinks)
MisLinks=rbind(c(56,51,3), MisLinks)
MisLinks=rbind(c(56,52,3), MisLinks)
MisLinks=rbind(c(56,53,3), MisLinks)
MisLinks=rbind(c(56,54,3), MisLinks)

MisLinks=rbind(c(72,49,3), MisLinks)
MisLinks=rbind(c(72,50,3), MisLinks)
MisLinks=rbind(c(72,51,3), MisLinks)
MisLinks=rbind(c(72,52,3), MisLinks)
MisLinks=rbind(c(72,53,3), MisLinks)
MisLinks=rbind(c(72,54,3), MisLinks)

MisNodes$group[56]=4

forceNetwork(Links = MisLinks, Nodes = MisNodes, 
             Source = "source", Target = "target", Value = "value", 
             NodeID = "name", Nodesize = "size", Group = "group",
             radiusCalculation = "Math.sqrt(d.nodesize)+6",
             opacity = 1.5, legend = T)
