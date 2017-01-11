exists_mis=0
exists_gap=0
fun.get_part_mis <- function(n,seq){
	#print(n)
   bit_add=1.2
  bit_h=3
  myh=h
  pv=0
	n<-lapply(n,as.numeric)[[1]]
	if(n>=1){
    for(i in 1:n){
    	#print("sssy")
      myh<-h
      pv<-0
    	line<-strsplit(readLines(con,n=1),"\t")[[1]]
      if(show_mis==1){
            l=lapply(line[1],as.numeric)[[1]]+1
            exists_mis<<-1
              if(l>1){
              for(i in 1:(l-1)){
                pv<-(pv+bit_add)
               if(pv>=96){pv<-0;myh<-(myh-bit_h)}
                   }
              }   
              #print(l)
               text(pv,myh,substr(seq,l,l),pos=4,cex=0.6)
      }
     } 
    }
}

fun.get_part_gap <- function(n){
	#print(n)
  bit_add=1.2
  bit_h=3
  myh=h
  pv=0
	n<-lapply(n,as.numeric)[[1]]
	if(n>=1 ){
    for(i in 1:n){
      myh<-h
      pv<-0
    	line<-strsplit(readLines(con,n=1),"\t")[[1]]
      if(show_gap==1){
                 l=lapply(line[1],as.numeric)[[1]]+1
                 exists_gap<<-1
              if(l>1){
              for(i in 1:(l-1)){
                pv<-(pv+bit_add)
               if(pv>=96){pv<-0;myh<-(myh-bit_h)}
                   }
              }   
               text(pv,myh,"|",pos=4,cex=0.6)
      }
    }
   }
}
fun.drawseq<- function(seq,s_e,color){
 	pv=0
 	l=lapply(s_e[1],as.numeric)[[1]]+1
 	r=lapply(s_e[2],as.numeric)[[1]]+l-1
  bit_add=1.2
  bit_h=3
 	myh=h
 	if(l>1){
 		for(i in 1:(l-1)){
      pv<-(pv+bit_add)
      if(pv>=96){pv<-0;myh<-(myh-bit_h)}
    }
 	}   
 	#print(seq)
    # print(count)
    for(i in l:r){
       #print(pv)
       if(pv>=96){
         myh<-(myh-bit_h)
       	 pv<-0
       }
   	   text(pv,myh,substr(seq,i,i),pos=4,col=color,cex=0.6)
       pv<-(pv+bit_add)
    }
}
fun.h_seq<- function(seq){
   
   t1<-nchar(seq)%/%80
   t2<-nchar(seq)%%80
   if(t2!=0)
    t1<-(t1+1)
    return(t1*3)
}
fun.h_seqdraw<- function(seq){
   l=nchar(seq)
   myh=h
   pv=0
   bit_add=1.2
   bit_h=3
   for(i in 1:l){
       #print(pv)
       if(pv>=96){
         myh<-(myh-bit_h)
         pv<-0
       }
       text(pv,myh,substr(seq,i,i),pos=4,col="grey",cex=0.6)
       pv<-(pv+bit_add)
    }
  
}
fun.get_part <- function(seq,color,n){
     for(i in 1:n){
         s_e<-strsplit(readLines(con,n=1),"\t")[[1]]
         fun.drawseq(seq,s_e,color)
         mis_n<-readLines(con,n=1) 
        # print(mis_n)
         fun.get_part_mis(mis_n,seq)
          gap_n<-readLines(con,n=1) 
         # print(gap_n)
         fun.get_part_gap(gap_n)
    }
}

seq_color1="grey"
seq_color2="grey"
fun.seqcolor <- function(col){
   if(col == "Blue") {seq_color1<<-"blue";seq_color2<<-"steelblue3"
   }
    if(col == "Green") {seq_color1<<-"green3";seq_color2<<-"palegreen3"
   }
    if(col == "Yellow") {seq_color1<<-"darkgoldenrod2";seq_color2<<-"lightgoldenrod2"
   }
   if(col == "Red") {seq_color1 <<-"#EE1A1A";seq_color2 <<-"HotPink2"
   }  
   if(col == "Purple") {seq_color1<<-"darkorchid3";seq_color2<<-"mediumpurple2"
   }    
   if(col == "Brown") {seq_color1<<-"chocolate3";seq_color2<<-"#C1B26A"
   }
    if(col == "Grey") { seq_color1<<-"gray14";seq_color2<<-"gray32"
   }
   if(col == "pink_blue") {seq_color1<<-"rosybrown3";seq_color2<<-"lightsteelblue3"
   }
    if(col == "red_gold") {seq_color1<<-"brown2";seq_color2<<-"goldenrod2"
   }
   if(col == "blue_green") {seq_color1<<-"#1F1F99";seq_color2<<-"#1F761F"
   }
   if(col == "green_yellow") {seq_color1<<-"#1F761F";seq_color2<<-"yellow3"
   }
   if(col == "red_green") {seq_color1<<-"#991F1F";seq_color2<<-"#1F761F"
   }
    if(col == "blue_yellow") {seq_color1<<-"#1F1F99";seq_color2<<-"yellow3"
   }
    if(col == "red_green") {seq_color1<<-"#991F1F";seq_color2<<-"#1F761F"
   }
}

Args <- commandArgs()
file_dir_seq=Args[6]
file_seq_format=Args[7]
file_output=Args[8]

pdf(file_output)
par(fig=c(0,1,0,1),font=1)
plot(1:100,1:100,type="n",xlab = "",ylab = "",axes = F)
#text(d_w,d_h,chrom_info[2],pos=d_p1)
format <- file(file_seq_format, "r")
n_line=1
line<-readLines(format,n=1)
show_mis=1
show_gap=1
while(length(line) != 0 ) {
     if(n_line==2){
      fun.seqcolor(line)
     }
    
      if(n_line==7){
        if(line=="no_mis"){show_mis<<-0}
     }
      if(n_line==8){
        if(line=="no_gap"){show_gap<<-0}
     }
     n_line<<-(n_line+1)
     line<-readLines(format,n=1)  

}
close(format)

#con <- file(file_dir_seq, "r")
con <- file(paste(file_dir_seq, "/info_seq.txt", sep=""), "r")
h=100
site_name<-readLines(con,n=1)
site_name<-paste("Circle Site: ", site_name,sep="")
text(0,h,site_name,pos=4,cex=0.8)
#two_c<-strsplit(readLines(con,n=1)," ")[[1]]
n<-readLines(con,n=1)
n<-lapply(n,as.numeric)[[1]]
h<-(h-10)
#print(n)
for(i in 1:n){
	
	name<-readLines(con,n=1)
	
	#print(name)
    strd<-readLines(con,n=1)
    if(strd=="+"){
       name<-paste(name,"  (strand: +)",sep="")
    }
   if(strd=="-"){
       name<-paste(name,"  (strand: -)",sep="")
    }
   # print(strd)
    seq<-readLines(con,n=1)
   
   # print(seq)
    part_num<<-strsplit(readLines(con,n=1),"\t")[[1]]
    #print(part_num)
    
    th<-fun.h_seq(seq)
    if((h-5-th)<0){
      par(fig=c(0,1,0,1),font=1)
      plot(1:100,1:100,type="n",xlab = "",ylab = "",axes = F)
      h<-100
    }
    text(0,h,name,pos=4,cex=0.7)
    h<<-(h-5)
    fun.h_seqdraw(seq) 
    fun.get_part(seq,seq_color1,lapply(part_num[1],as.numeric)[[1]])
    fun.get_part(seq,seq_color2,lapply(part_num[2],as.numeric)[[1]])
    h<<-(h-th-2)
}
close(con)
if(exists_mis==1){
  text(0,h,"black color: mismatch",pos=4,cex=0.7)
   h<<-(h-3)
}
if(exists_gap==1){
  text(0,h,"|: gap",pos=4,cex=0.7)
}
#text(30,30,S[2],col="black",pos=4)
#text(30+2.6315*15,30,"BCC",col="red",pos=4)

