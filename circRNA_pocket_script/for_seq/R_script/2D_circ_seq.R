
b_firstp=0
my_size=0.8

Args <- commandArgs()

file_dir=Args[6]
file_format=Args[7]
file_output=Args[8]

con_gene <- file(paste(file_dir, "\\info.txt", sep=""), "r")
con_format <- file(file_format, "r")

con_gene2 <- file(paste(file_dir, "\\info.txt", sep=""), "r")
line<-readLines(con_gene2,n=1)
line<-readLines(con_gene2,n=1)
line<-readLines(con_gene2,n=1)
site_arr<-strsplit(line,"_")[[1]]
site_long<-(lapply(site_arr[3],as.numeric)[[1]]-lapply(site_arr[2],as.numeric)[[1]])
close(con_gene2)

fun.seqdraw <- function(c,s,p,w){
  
   ss=array(0,dim=c(3,1))
    ss[1]=lapply(s[1],as.numeric)[[1]]
    ss[2]=lapply(s[2],as.numeric)[[1]]
    ss[3]=lapply(s[3],as.numeric)[[1]]

  
    if(b_firstp==0){par(fig=c(0,0.9,0.05,0.95));b_firstp<<-1}
    else{par(new=TRUE,fig=c(0,0.9,0.05,0.95))}
    mcol=array("",dim=c(3,0))
    mcol[2]=c[1]
    pie(ss, col= mcol,radius = p,xlab = "",ylab = "",border=0,labels = "",clockwise=TRUE, init.angle = 90 )
  
   par(new=TRUE,fig=c(0,0.9,0.05,0.95))
   pie(1, col = "white",radius = p-w,xlab = "",ylab = "",border=0,labels = "")
}

fun.seqline <- function(c,s){
   ss=array(0,dim=c(2,1))
    ss[1]=lapply(s[1],as.numeric)[[1]]
    ss[2]=lapply(s[2],as.numeric)[[1]]
    
     #par(new=TRUE,fig=c(0,1,0,1),font=1)
  seqpos1=ss[1]*seq_line_long/site_long
  seqpos2=(ss[1]+ss[2])*seq_line_long/site_long
 
    rect(seq_line_posw+seqpos1, seq_line_posh, seq_line_posw+seqpos2, seq_line_posh+seq_line_wid,col=c[1],border=0)
   
}

fun.seqdirec <- function(s,seql){
   pos=seq_line_posw+seq_line_long*seql/site_long
   long=2
    x <- seq(pos+long,pos+long+2,length = 1000);
    y <- seq(seq_line_posh,seq_line_posh,length = 1000);
    lines(x, y)
    if(s=="+"){
         x <- seq(pos+long+1,pos+long+2,length = 1000);
         y <- seq(seq_line_posh+0.8,seq_line_posh,length = 1000);
         lines(x, y)
    }
    if(s=="-"){
         x <- seq(pos+long+1,pos+2,length = 1000);
         y <- seq(seq_line_posh+0.8,seq_line_posh,length = 1000);
         lines(x, y)
    }

}
line_wid=0.01
#mis_color=array("grey10",dim=c(2,1))
mis_color=array("black",dim=c(1,1))
#gap_color=array("grey90",dim=c(2,1))
gap_color=array("white",dim=c(1,1))
seq_color1=array("grey",dim=c(1,1))
seq_color2=array("grey",dim=c(1,1))
seq_pos=0.8
seq_wid=0.05
seq_interv=0.05
show_mis=1
show_gap=1
line=0
seq_line_long=50
seq_line_wid=2
seq_line_posw=0
seq_line_posh=20
seq_line_interv=6
fun.seq <- function(seq_c,num){  
      line<-readLines(con_seq,n=1) 
      match_struct=strsplit(line,"\t")[[1]]
      if(num==1){ fun.seqdraw(seq_c,match_struct,seq_pos,seq_wid)}
     else{fun.seqline(seq_c,match_struct)}
      #fun.seqline(seq_c,match_struct,posw)
      line<-readLines(con_seq,n=1) 
      n_mis=lapply(line,as.numeric)[[1]]
      if(n_mis>0){
          for(k in 1:n_mis){
           line<-readLines(con_seq,n=1) 
           if(show_mis==1){
                  mis_struct=strsplit(line,"\t")[[1]]
                  if(num==1){fun.seqdraw(mis_color,mis_struct,seq_pos,seq_wid)}
                  else{fun.seqline(mis_color,mis_struct)}
                 # fun.seqline(mis_color,mis_struct,posw)
            }        
          }
      }
    
        line<-readLines(con_seq,n=1) 
        n_gap=lapply(line,as.numeric)[[1]]
      if(n_gap>0){
           for(k in 1:n_gap){
           line<-readLines(con_seq,n=1) 
           if(show_gap==1){
                gap_struct=strsplit(line,"\t")[[1]]
                if(num==1){fun.seqdraw(gap_color,gap_struct,seq_pos,seq_wid)}
                else{fun.seqline(gap_color,gap_struct)}
                # fun.seqline(gap_color,gap_struct,posw)
           }         
        }
      }    
}
fun.seqcolor <- function(col){
    if(col == "Blue") {seq_color1[1]<<-"blue";seq_color2[1]<<-"steelblue3"
   }
    if(col == "Green") {seq_color1[1]<<-"green3";seq_color2[1]<<-"palegreen3"
   }
    if(col == "Yellow") {seq_color1[1]<<-"darkgoldenrod2";seq_color2[1]<<-"lightgoldenrod2"
   }
   if(col == "Red") {seq_color1[1] <<-"#EE1A1A";seq_color2[1] <<-"HotPink2"
   }  
   if(col == "Purple") {seq_color1[1]<<-"darkorchid3";seq_color2[1]<<-"mediumpurple2"
   }    
   if(col == "Brown") {seq_color1[1]<<-"chocolate3";seq_color2[1]<<-"#C1B26A"
   }
    if(col == "Grey") { seq_color1[1]<<-"gray14";seq_color2[1]<<-"gray32"
   }
   if(col == "pink_blue") {seq_color1[1]<<-"rosybrown3";seq_color2[1]<<-"lightsteelblue3"
   }
    if(col == "red_gold") {seq_color1[1]<<-"brown2";seq_color2[1]<<-"goldenrod2"
   }
   if(col == "blue_green") {seq_color1[1]<<-"#1F1F99";seq_color2[1]<<-"#1F761F"
   }
   if(col == "green_yellow") {seq_color1[1]<<-"#1F761F";seq_color2[1]<<-"yellow3"
   }
   if(col == "red_green") {seq_color1[1]<<-"#991F1F";seq_color2[1]<<-"#1F761F"
   }
    if(col == "blue_yellow") {seq_color1[1]<<-"#1F1F99";seq_color2[1]<<-"yellow3"
   }
    if(col == "blue_red") {seq_color1[1]<<-"#1F1F99";seq_color2[1]<<-"#991F1F"
   }

}

exon_b=0;intron_b=0;inter_b=0;UTR_b=0

exon_color="blue"
intron_color="cyan2"
UTR_color="steelblue3"
inter_color="grey40"
shadow_color="grey80"


 be_link=0.1
fun.color <- function(col){
     if(col == "Red") {
      exon_color <<-"#EE1A1A";intron_color <<-"lightpink2";UTR_color <<-"HotPink2";
   }
   if(col == "Yellow") {
      exon_color<<-"darkgoldenrod2";intron_color<<-"#F2F200";UTR_color<<-"lightgoldenrod2";
   }
   if(col == "Green") {
      exon_color<<-"green3";intron_color<<-"darkolivegreen2";UTR_color<<-"palegreen3";
   }
   if(col == "Purple") {
      exon_color<<-"darkorchid3";intron_color<<-"#EF84EA";UTR_color<<-"mediumpurple2";
   }    
   if(col == "Brown") {
      exon_color<<-"chocolate3";intron_color<<-"#F09B49";UTR_color<<-"#C1B26A";#darkkhaki
   }
   if(col == "Grey") {
      exon_color<<-"gray14";intron_color<<-"gray50";UTR_color<<-"gray32";
   }
   if(col == "pink_white_blue") {
      exon_color<<-"rosybrown3";intron_color<<-"antiquewhite3";UTR_color<<-"lightsteelblue3";
   }
    if(col == "red_wood_gold") {
      exon_color<<-"brown2";intron_color<<-"burlywood2";UTR_color<<-"goldenrod2";
    }
   if(col == "blue_red_green") {
      exon_color<<-"darkblue";intron_color<<-"darkred";UTR_color<<-"darkgreen";
   }
    if(col == "blue_yellow_green") {
      exon_color<<-"navy";intron_color<<-"yellow3";UTR_color<<-"#2A8F2A";#forestgreen
   }
}

UTR_flag=1

inter_r=0.95
inter_l=0.85
exon_r=1
exon_l=0.8
narrow_tag=1
line_cut="S"
fun.width <- function(wid){
  line_cut<<-strsplit(wid,"_")[[1]]
  my_width=lapply(line_cut[1],as.numeric)[[1]]
  my_interval=lapply(line_cut[2],as.numeric)[[1]]
      inter_r<<-(1-my_width/200)
      inter_l<<-(1-my_width/200*3)
      exon_r<<-1
      exon_l<<-(1-my_width/50)
      be_link<<-((my_width+my_interval)/50)
}

my_clockwise = FALSE
d_w=82;d_h=50;d_p1=3;d_p2=1
x11=75;x12=89;y11=50;y12=50;
fun.angle <- function(){
   #if(clok == "clockwise") {
      d_w<<-44;d_h<<-50;d_p1<<-4;d_p2<<-4
      x11<<-44;x12<<-44;y11<<-50;y12<<-92;      
       my_clockwise<<-TRUE
   #}  
}

fun.direction<- function(){
  if(gene_direction!="+" && gene_direction!="-"){return}
   draw_directon=TRUE
   if(my_clockwise==FALSE){
     draw_directon<-(!draw_directon)
   }
   if(gene_direction=="-"){
    draw_directon<-(!draw_directon)
   }
   pos=(exon_r-be_link)*my_size
     s=array(50,dim=c(3,1)) 
     s[2]=20 
      c=array("white",dim=c(3,1)) 
      c[2]="grey10"
    par(new=TRUE,fig=c(0,0.9,0.05,0.95))
     pie(s, col = c,radius = pos,xlab = "",ylab = "",border=0,labels = "",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)
     par(new=TRUE,fig=c(0,0.9,0.05,0.95))
     pie(1, col = "white",radius = pos-0.008,xlab = "",ylab = "",border=0,labels = "",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)
    
     par(new=TRUE,fig=c(0,1,0,1),font=1)
     plot(1:100,1:100,type="n",xlab = "",ylab = "",axes = F)
    #direction_pos=(exon_r-be_link)*my_size
    if(draw_directon){
      temp<-5
      x <- seq(44-0.5*pos*45,44-0.5*pos*45+0.8660254*temp,length = 1000);
      y <- seq(51-0.8660254*pos*45,51-0.8660254*pos*45+0.5*temp,length = 1000);
       lines(x, y)
    }
     else{
       temp<-5
      x <- seq(44+0.5*pos*45*0.9,44+0.5*pos*45*0.9-0.8660254*temp,length = 1000);
      y <- seq(51-0.8660254*pos*45*0.95,51-0.8660254*pos*45*0.95+0.5*temp,length = 1000);
       lines(x, y)
    }

}
#pdf("cir.pdf")
#pdf("mymy.pdf")
pdf(file_output)


n_line=0
count_circ=0
ios_name=""
ios_name_line=""
chrom_info=""
chrom_info_line=""
gene_direction="0"
show_gene_direction=FALSE
line<-readLines(con_gene,n=1)
while(length(line) != 0 ) {
     if(n_line==0){    
      count_circ<<- lapply(line,as.numeric)[[1]]
      #fun.color(line)
     }
     if(n_line==1){
      #ios_name_line<<-line
      ios_name<<-strsplit(line,"\t")[[1]]
       if(length(ios_name)>2){
         ios_name_line<<-paste(ios_name[1], " ",ios_name[2],"...",sep="")
       }
       else{
          for(i in 1:length(ios_name)){
            ios_name_line<<-paste(ios_name_line, " ",ios_name[i],sep="")
          }
       }
         if(nchar(ios_name_line)>50){
           ios_name_line<<-ios_name[1]
           if(length(ios_name)>2){ ios_name_line<<-paste(ios_name[1], "...",sep="")}
           if(nchar(ios_name_line)>50){ios_name_line<<-"..."}
         }
     }
     if(n_line==2){
      chrom_info_line<<-line
      chrom_info<<-strsplit(line,"_")[[1]]
     }
     if(n_line==3){
      gene_direction<<-line
     }
     n_line<<-(n_line+1)
     line<-readLines(con_gene,n=1)    
}
close(con_gene)

n_line<-1
line<-readLines(con_format,n=1)
fun.angle()
while(length(line) != 0 ) {

     if(n_line==1){
      fun.color(line)
     }
      if(n_line==2){
      fun.seqcolor(line)
     }
     if(n_line==3){
      circle_num_type<-strsplit(line,"\t")[[1]]
      if(as.numeric(count_circ)<2){fun.width(circle_num_type[1])}
      else{fun.width(circle_num_type[2])}
      #fun.width(line)
     }
      if(n_line==4){
      seq_pw<<-strsplit(line,"_")[[1]]
      seq_wid<<-lapply(seq_pw[1],as.numeric)[[1]]/50
      seq_interv<<-lapply(seq_pw[2],as.numeric)[[1]]/50
     }

     if(n_line==5){
      if(line=="no_UTR"){UTR_flag<<-0}
     }
     if(n_line==6){
       if(line=="show_gene_direction"){
          show_gene_direction<<-TRUE
       }
     }
      if(n_line==7){
        if(line=="no_mis"){show_mis<<-0}
     }
      if(n_line==8){
        if(line=="no_gap"){show_gap<<-0}
     }
      if(n_line==9){
        my_size<<-lapply(line,as.numeric)[[1]]
     }
     n_line<<-(n_line+1)
     line<-readLines(con_format,n=1)   
}
close(con_format)


#con_seq <- file("F:\\d\\circRNA_pocket\\Release\\output\\R_result\\Chr1_378967_379205\\info_seqc.txt", "r")
con_seq <- file(paste(file_dir, "\\info_seqc.txt", sep=""), "r")
#A=array(0,dim=c(1000,1))
#line<-readLines(con,n=1) 


#plot(1:100,1:100,type="n",xlab = "",ylab = "")
#par(new=TRUE,fig=c(0,1,0,1),font=1)
seq_n<-lapply(readLines(con_seq,n=1) ,as.numeric)[[1]]
seq_pos=my_size
seq_pos<-(seq_pos+(seq_wid+seq_interv)*(seq_n+1))
for(i in 1:seq_n){
   #seq_line_posh<<-(seq_line_posh+seq_line_wid+seq_line_interv)
   seq_pos<<-(seq_pos-(seq_wid+seq_interv))
   line<-readLines(con_seq,n=1) 
   line<-readLines(con_seq,n=1)
   #text(seq_line_posw,seq_line_posh+seq_line_wid+1,line, cex = 0.8,pos=3)
   line<-readLines(con_seq,n=1) 
   ios_num<<-strsplit(line,"\t")[[1]]
   f_s<-lapply(ios_num[1],as.numeric)[[1]]
   for(j in 1:f_s){
          fun.seq(seq_color1,1)     
   }
    e_s<-lapply(ios_num[2],as.numeric)[[1]]
    for(j in 1:e_s){
          fun.seq(seq_color2,1)     
   }

}
close(con_seq)

















#fun.test(a = 10, b = 8, method = "add")

#Args <- commandArgs()

#file_dir1=Args[6]
#file_dir2=Args[7]
#file_dir3=Args[8]

#file_dir1="F:\\d\\circRNA_pocket\\Release\\output\\R_result\\Chr5_21715716_21724565"
#file_dir2="F:\\d\\circRNA_pocket\\Release\\R_script\\Setting\\circ_format.txt"
#file_dir3="itis.pdf"
#con <- file(paste(file_dir1, "\\info.txt", sep=""), "r")
#if(is.na(Args[9])){  
#   con <- file(paste(file_dir1, "\\info.txt", sep=""), "r")
#  }else{
#    con <- file(paste(file_dir1, "\\",Args[9], sep="") , "r")
#  }

#con2 <- file(file_dir2, "r")





for(p in 1:as.numeric(count_circ)){

browers <- read.table(paste(file_dir, "\\",ios_name[p], sep="") ,stringsAsFactors = FALSE)
ty=browers[,1]

if(p>1){
  inter_r<<-(inter_r-be_link)
    inter_l<<-(inter_l-be_link)
    exon_r<<-(exon_r-be_link)
    exon_l<<-(exon_l-be_link)
}

for(j in 1:length(ty)){

    if(ty[j]=="CDS"){
      ty[j]=exon_color      
      exon_b=1
    }
    if(ty[j] == "intron"){
      ty[j]=intron_color
      intron_b=1
    }
    
    if(ty[j] == "intergenic"){
      ty[j]="0"
      inter_b=1
    }
    if(ty[j]=="five_prime_UTR" || ty[j]=="three_prime_UTR" || ty[j]=="UTR"){
      if(UTR_flag==1){
          ty[j]=UTR_color
          UTR_b=1
      }
      else{
          ty[j]=exon_color      
           exon_b=1
      }
    }
}

y=browers[,2]
x=1

#if((show_gene_direction &&(gene_direction=="+" || gene_direction=="-"))|| p>1){
  par(new=TRUE,fig=c(0,0.9,0.05,0.95))
#}
pie(y, col = inter_color,radius = inter_r*my_size,xlab = "",ylab = "",border=0, labels = "")
par(new=TRUE,fig=c(0,0.9,0.05,0.95))
pie(y, col = "white",radius = inter_l*my_size,xlab = "",ylab = "",border=0,labels = "")

#par(new=TRUE,cex.main=1)
par(new=TRUE,fig=c(0,0.9,0.05,0.95),font.main=4)
if (p==1){
  pie(y, col = ty,density = NULL, lty = NULL, angle = 45,main = paste(chrom_info[1], ":",ios_name_line, sep="")  ,radius = exon_r*my_size,xlab = "",ylab = "",border=0, labels = "",col.main="black",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0,)
}
else{
  pie(y, col = ty, angle = 45,radius = exon_r*my_size,xlab = "",ylab = "",border=0,labels = "",density = NULL,lty = NULL,clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)
}
par(new=TRUE,fig=c(0,0.9,0.05,0.95))
pie(x, col = "white",radius = exon_l*my_size,xlab = "",ylab = "",border=0,labels = "")

}

if(show_gene_direction){
  fun.direction()
}

par(new=TRUE,fig=c(0,1,0,1),font=1)
plot(1:100,1:100,type="n",xlab = "",ylab = "",axes = F)
#plot(1:100,1:100,type="n",xlab = "",ylab = "")

gene_line_long=25
gene_line_wid=2
gene_line_posw=60
gene_line_posh=20
gene_line_interv=6



add=(exon_l-(exon_r-exon_l)*0.05*narrow_tag)*my_size*45
text(d_w,d_h+add-2,chrom_info[2],pos=d_p1,srt=270)
text(d_w-4.5,d_h+add-2,chrom_info[3],pos=d_p2,srt=270)

#x <- seq(x11,x12,length = 1000);
#y <- seq(y11,y12,length = 1000);
x <- seq(44,44,length = 1000);
y <- seq(50,d_h+add,length = 1000);
lines(x, y)

mi_h=3;mi_m=1.2;sta_p=30;sta_l=80+5;m_l=4;m_shadow=0.2
if (exon_b==1){
  rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =exon_color,border=0,)
  if(UTR_flag==1){
    text(sta_l+m_l,sta_p-mi_h*0.6,"CDS", cex = 0.8,pos=4)
  }
  else{
    text(sta_l+m_l,sta_p-mi_h*0.6,"Exon", cex = 0.8,pos=4)
  }
  sta_p<<-(sta_p-mi_h-mi_m)
}
if (intron_b==1){
  rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =intron_color,border=0,)
  text(sta_l+m_l,sta_p-mi_h*0.6,"Intron", cex = 0.8,pos=4)
  sta_p<<-(sta_p-mi_h-mi_m)
}
if (inter_b==1){
  rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =inter_color,border=0,)
  text(sta_l+m_l,sta_p-mi_h*0.6,"Intergenic", cex = 0.8,pos=4)
  sta_p<<-(sta_p-mi_h-mi_m)
}
if (UTR_b==1){
  rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =UTR_color,border=0,)
  text(sta_l+m_l,sta_p-mi_h*0.6,"UTR", cex = 0.8,pos=4)
  sta_p<<-(sta_p-mi_h-mi_m)
}
if(show_gene_direction && (gene_direction=="-" || gene_direction=="+")){
   x <- seq(sta_l,sta_l+m_l,length = 1000);
  y <- seq(sta_p-mi_h*0.8,sta_p-mi_h*0.8,length = 1000);
  lines(x, y)
 
  x <- seq(sta_l+m_l*0.7,sta_l+m_l,length = 1000);
  y <- seq(sta_p-mi_h*0.5,sta_p-mi_h*0.8,length = 1000);
  lines(x, y)
 
   text(sta_l+m_l,sta_p-mi_h*0.6,"Strand", cex = 0.8,pos=4)
    sta_p=sta_p-mi_h-mi_m
}

par(fig=c(0,1,0,1),font=1)
plot(1:100,1:100,type="n",xlab = "",ylab = "",axes = F)
#con_seq <- file("F:\\d\\circRNA_pocket\\Release\\output\\R_result\\Chr1_378967_379205\\info_seql.txt", "r")
con_seq <- file(paste(file_dir, "\\info_seql.txt", sep=""), "r")

seq_n<-lapply(readLines(con_seq,n=1) ,as.numeric)[[1]]
seq_line_posh<<-(seq_line_posh+(seq_line_wid+seq_line_interv)*(seq_n+1))
for(i in 1:seq_n){
   seq_line_posh<<-(seq_line_posh-(seq_line_wid+seq_line_interv))
   #seq_pos<<-(seq_pos+seq_wid+seq_interv)
   line<-readLines(con_seq,n=1) 
   text(seq_line_posw-2,seq_line_posh+seq_line_wid+2.5,line, cex = 0.8,pos=4)
   line<-readLines(con_seq,n=1) 
   seql=lapply(line,as.numeric)[[1]]
   line<-readLines(con_seq,n=1) 
   fun.seqdirec(line,seql)  
   rect(seq_line_posw, seq_line_posh, seq_line_posw+seq_line_long*seql/site_long, seq_line_posh+seq_line_wid,col ="grey",border=0)
   line<-readLines(con_seq,n=1) 
   ios_num<<-strsplit(line,"\t")[[1]]
   f_s<-lapply(ios_num[1],as.numeric)[[1]]
   for(j in 1:f_s){
          fun.seq(seq_color1,2)     
   }
    e_s<-lapply(ios_num[2],as.numeric)[[1]]
    for(j in 1:e_s){
          fun.seq(seq_color2,2)     
   }
}
close(con_seq)

if(show_gene_direction && (gene_direction=="-" || gene_direction=="+")){
    x <- seq(gene_line_posw+0.3*gene_line_long,gene_line_posw,length = 1000);
    y <- seq(17,17,length = 1000);
    lines(x, y)
    if(gene_direction=="+"){
         x <- seq(gene_line_posw+0.2*gene_line_long,gene_line_posw+0.3*gene_line_long,length = 1000);
         y <- seq(18,17,length = 1000);
         lines(x, y)
    }
    if(gene_direction=="-"){
         x <- seq(gene_line_posw+0.1*gene_line_long,gene_line_posw,length = 1000);
         y <- seq(18,17,length = 1000);
         lines(x, y)
    }
}




fun.geneline <- function(y,ty,gt,gb){
  l=length(y)
   seql=0
  for(j in 1:l){seql<-(seql+y[j])}
  rate=gene_line_long/seql
  pre=0;temp=0
   for(j in 1:l){
     temp<-(temp+y[j])
     rect(gene_line_posw+pre*rate, gt, gene_line_posw+temp*rate, gb,col=ty[j],border=0)
     pre<-(pre+y[j])
  }

}


gene_line_cex=0.8
if(as.numeric(count_circ)>9){
  gene_line_interv<-4.5
  gene_line_wid<-1.8
  gene_line_cex<-0.7
}

for(p in as.numeric(count_circ):1){

browers <- read.table(paste(file_dir, "\\",ios_name[p], sep="") ,stringsAsFactors = FALSE)

gene_line_posh<<-(gene_line_posh+gene_line_wid+gene_line_interv)
text(gene_line_posw-2,gene_line_posh+gene_line_wid+2,ios_name[p], cex = gene_line_cex,pos=4)

ty=browers[,1]
for(j in 1:length(ty)){

    if(ty[j]=="CDS"){
      ty[j]=exon_color      
      exon_b=1
    }
    if(ty[j] == "intron"){
      ty[j]=intron_color
      intron_b=1
    }
    
    if(ty[j] == "intergenic"){
      ty[j]="0"
      inter_b=1
    }
    if(ty[j]=="five_prime_UTR" || ty[j]=="three_prime_UTR" || ty[j]=="UTR"){
      if(UTR_flag==1){
          ty[j]=UTR_color
          UTR_b=1
      }
      else{
          ty[j]=exon_color      
           exon_b=1
      }
    }
}

y=browers[,2]
x=1

rect(gene_line_posw, gene_line_posh+gene_line_wid*0.15, gene_line_posw+gene_line_long,gene_line_posh+gene_line_wid*0.85,col =inter_color,border=0,)
fun.geneline(y,ty,gene_line_posh,gene_line_posh+gene_line_wid)

}


text(5,10,"sequence", cex = 0.8,pos=4)
text(gene_line_posw,10,"gene struct", cex = 0.8,pos=4)
text(0,0,"(up to down,outer to inner)", cex = 0.7,pos=4)





