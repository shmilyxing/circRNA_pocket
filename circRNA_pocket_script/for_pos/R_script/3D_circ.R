exon_b=0;intron_b=0;inter_b=0;UTR_b=0

exon_color="blue"
intron_color="cyan2"
UTR_color="steelblue3"
inter_color="grey40"
shadow_color="grey80"
exon_color_low="#2E45F6"
intron_color_low="cyan"
UTR_color_low="#58A6E7"
inter_color_low="grey48"
exon_color_high="#0000A7"
intron_color_high="cyan3"
UTR_color_high="#3F76A4"
inter_color_high="grey36"

 be_link=0.1
fun.color <- function(col){
   if(col == "Red") {
      exon_color <<-"#EE1A1A";exon_color_low <<-"#FF1313";exon_color_high <<-"#CD0F0F"
      intron_color <<-"lightpink2";intron_color_low <<-"lightpink1";intron_color_high <<-"lightpink3"
      UTR_color <<-"HotPink2";UTR_color_low <<-"HotPink1";UTR_color_high <<-"HotPink3"
   }
   if(col == "Yellow") {
      exon_color<<-"darkgoldenrod2";exon_color_low<<-"darkgoldenrod1";exon_color_high<<-"darkgoldenrod3"
      intron_color<<-"#F2F200";intron_color_low<<-"yellow";intron_color_high<<-"#D8D800"
      UTR_color<<-"lightgoldenrod2";UTR_color_low<<-"#F8E687";UTR_color_high<<-"lightgoldenrod3"
   }
   if(col == "Green") {
      exon_color<<-"green3";exon_color_low<<-"#00DD00";exon_color_high<<-"#00A400"
      intron_color<<-"darkolivegreen2";intron_color_low<<-"darkolivegreen1";intron_color_high<<-"darkolivegreen3"
      UTR_color<<-"palegreen3";UTR_color_low<<-"palegreen2";UTR_color_high<<-"palegreen4"
   }
   if(col == "Purple") {
      exon_color<<-"darkorchid3";exon_color_low<<-"#B039EB";exon_color_high<<-"darkorchid4"
      intron_color<<-"#EF84EA";intron_color_low<<-"#FF8DFA";intron_color_high<<-"#D175CD"
      UTR_color<<-"mediumpurple2";UTR_color_low<<-"mediumpurple1";UTR_color_high<<-"mediumpurple3"
   }    
   if(col == "Brown") {
      exon_color<<-"chocolate3";exon_color_low<<-"#DF6E1F";exon_color_high<<-"#A55217"
      intron_color<<-"#F09B49";intron_color_low<<-"tan1";intron_color_high<<-"tan3"
     UTR_color<<-"#C1B26A";UTR_color_low<<-"#C9BC7E";UTR_color_high<<-"#A29453"#darkkhaki
   }
    if(col == "Grey") {
      exon_color<<-"gray14";exon_color_low<<-"gray18";exon_color_high<<-"gray8"
      intron_color<<-"gray50";intron_color_low<<-"gray55";intron_color_high<<-"gray40"
      UTR_color<<-"gray32";UTR_color_low<<-"gray38";UTR_color_high<<-"gray28"
   }
   
    if(col == "pink_white_blue") {
      exon_color<<-"rosybrown3";exon_color_low<<-"#E5ADAD";exon_color_high<<-"#AF8484"
      intron_color<<-"antiquewhite3";intron_color_low<<-"#E5D6C4";intron_color_high<<-"#AFA496"
      UTR_color<<-"lightsteelblue3";UTR_color_low<<-"#B5CAE5";UTR_color_high<<-"#8A9AAF"
   }
    if(col == "red_wood_gold") {
      exon_color<<-"brown2";exon_color_low<<-"brown1";exon_color_high<<-"brown3"
      intron_color<<-"burlywood2";intron_color_low<<-"burlywood1";intron_color_high<<-"burlywood3"
      UTR_color<<-"goldenrod2";UTR_color_low<<-"goldenrod1";UTR_color_high<<-"goldenrod3"
   }
   if(col == "blue_red_green") {
      exon_color<<-"darkblue";exon_color_low<<-"#1F1F99";exon_color_high<<-"#000052"
      intron_color<<-"darkred";intron_color_low<<-"#991F1F";intron_color_high<<-"#520000"
      UTR_color<<-"darkgreen";UTR_color_low<<-"#1F761F";UTR_color_high<<-"#003B00"
   }
    if(col == "blue_yellow_green") {
      exon_color<<-"navy";exon_color_low<<-"#1F1F8F";exon_color_high<<-"#00004B"
      intron_color<<-"yellow3";intron_color_low<<-"#D3D31F";intron_color_high<<-"#B3B300"
      UTR_color<<-"#2A8F2A";UTR_color_low<<-"#3D993D";UTR_color_high<<-"#186618"#forestgreen
   }
}

inter_r=0.95
inter_l=0.85
exon_r=1
exon_l=0.8
my_size=0.8
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
fun.angle <- function(clok){
   if(clok == "clockwise") {
     	d_w<<-44;d_h<<-88;d_p1<<-4;d_p2<<-2
      x11<<-44;x12<<-44;y11<<-85;y12<<-92;      
       my_clockwise<<-TRUE
   }  
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
   pos=0.9
     s=array(50,dim=c(3,1)) 
     s[2]=20 
      c=array("white",dim=c(3,1)) 
      c[2]="grey10"
     par(fig=c(0,0.9,0.05,0.95))
     pie(s, col = c,radius = pos,xlab = "",ylab = "",border=0,labels = "",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)
     par(new=TRUE,fig=c(0,0.9,0.05,0.95))
     pie(1, col = "white",radius = pos-0.005,xlab = "",ylab = "",border=0,labels = "",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)
    
     par(new=TRUE,fig=c(0,1,0,1),font=1)
     plot(1:100,1:100,type="n",xlab = "",ylab = "",axes = F)
     #plot(1:100,1:100,type="n",xlab = "",ylab = "")
     r=39
    if(draw_directon){
       #  x <- seq(pic_mid-0.5*direction_pos*45,pic_mid-0.5773503*direction_pos*45,length = 100)
     #   y <- seq(pic_mid,pic_mid,length = 100)
       if(my_clockwise){
             x <- seq(44-0.5*pos*r,44-0.5*pos*r+2.8,length = 1000);
             y <- seq(50-0.8660254*pos*r-2,50-0.8660254*pos*r-5.8,length = 1000);           
        }
        else{
            x <- seq(44-0.8660254*pos*r+0.6,44-0.8660254*pos*r-2.7,length = 1000);   
             y <- seq(50+0.5*pos*r+2,50+0.5*pos*r-1.1,length = 1000);  
        }
    
     }
     else{
        if(my_clockwise){
             x <- seq(44+0.5*pos*r,44+0.5*pos*r-2.8,length = 1000);
             y <- seq(50-0.8660254*pos*r-2,50-0.8660254*pos*r-5.8,length = 1000);
        }
        else{
             x <- seq(44-0.8660254*pos*r+0.6,44-0.8660254*pos*r-3,length = 1000);   
             y <- seq(50-0.5*pos*r-1.3,50-0.5*pos*r+1.8,length = 1000);    
       }
     }
      lines(x, y)
}

UTR_flag=1
#fun.test(a = 10, b = 8, method = "add")

Args <- commandArgs()
#file_dir=Args[6]
file_dir1=Args[6]
file_dir2=Args[7]
file_dir3=Args[8]
con=""
if(is.na(Args[9])){  
    con <- file(paste(file_dir1, "/info.txt", sep=""), "r")
  }else{
    con <- file(paste(file_dir1, "/",Args[9], sep="") , "r")
  }
#file_dir1="F:\\d\\circRNA_pocket\\Release\\output\\R_result\\Chr5_21715716_21724565"
#file_dir2="F:\\d\\circRNA_pocket\\Release\\R_script\\Setting\\circ_format.txt"
#file_dir3="itis.pdf"


con2 <- file(file_dir2, "r")

n_line=0
count_circ=0
ios_name=""
ios_name_line=""
chrom_info=""
chrom_info_line=""
gene_direction="0"
line=readLines(con,n=1)
while(length(line) != 0 ) {
     if(n_line==0){
     	count_circ<<-line
     	fun.color(line)
     }
     if(n_line==1){
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
     line<-readLines(con,n=1)    
}
close(con)
n_line<<-0
line=readLines(con2,n=1)
show_gene_direction=FALSE
while(length(line) != 0 ) {

     if(n_line==0){
     	fun.color(line)
     }
     if(n_line==1){
     	circle_num_type<-strsplit(line,"\t")[[1]]
     	if(as.numeric(count_circ)<2){fun.width(circle_num_type[1])}
     	else{fun.width(circle_num_type[2])}
     	#fun.width(line)
     }
     if(n_line==2){fun.angle(line)}
     if(n_line==3){
     	if(line=="no_UTR"){UTR_flag<<-0}
     }
      if(n_line==4){
       if(line=="show_gene_direction"){
          show_gene_direction<<-TRUE
       }
     }
     n_line<<-(n_line+1)
     line<-readLines(con2,n=1)   
}
close(con2)


#pdf(paste(file_dir, "\\", chrom_info_line,".pdf",sep=""))
pdf(file_dir3)

ramp_exon <- colorRamp(c(exon_color, exon_color_high)) 
ramp_intron <- colorRamp(c(intron_color, intron_color_high)) 
ramp_UTR <- colorRamp(c(UTR_color, UTR_color_high))
ramp_inter<-colorRamp(c(inter_color, inter_color_high)) 
my_size=0.8
par(fig=c(0,0.9,0.05,0.95))

if(show_gene_direction){
  fun.direction()
}
for(p in 1:as.numeric(count_circ)){

browers <- read.table(paste(file_dir1, "/",ios_name[p], sep="") ,stringsAsFactors = FALSE)

if(p>1){
	inter_r<<-(inter_r- be_link)
    inter_l<<-(inter_l- be_link)
    exon_r<<-(exon_r- be_link)
    exon_l<<-(exon_l- be_link)
}

ty_name=browers[,1]
ty=browers[,1]
ty_low=browers[,1]

for(j in 1:length(ty)){
	  if(ty[j]=="CDS"){
			ty[j]=exon_color
			ty_low[j]=exon_color_low
			exon_b=1
		}
#"#2E45F6"
		if(ty[j] == "intron"){
			ty[j]=intron_color
			ty_low[j]=intron_color_low
			intron_b=1
		}
		if(ty[j]=="five_prime_UTR" || ty[j]=="three_prime_UTR" || ty[j]=="UTR"){
			if(UTR_flag==1){
				ty[j]=UTR_color
			    ty_low[j]=UTR_color_low
			    UTR_b=1
			}
			else{
				ty[j]=exon_color
			    ty_low[j]=exon_color_low
			    exon_b=1
			}
		
		}
		if(ty[j] == "intergenic"){
			ty[j]="0"
			ty_low[j]="0"
			#cat(rainbow(n))
			inter_b=1
		}		
}


y=browers[,2]
x=1
if((show_gene_direction &&(gene_direction=="+" || gene_direction=="-"))|| p>1){
  par(new=TRUE,fig=c(0,0.9,0.05,0.95))
}
pie(1, col = inter_color_low,radius = inter_r*my_size,xlab = "",ylab = "",border=0, labels = "")
par(new=TRUE,fig=c(0,0.9,0.05,0.95))
pie(1, col = inter_color,radius = (inter_r-(inter_r-inter_l)*0.2)*my_size,xlab = "",ylab = "",border=0, labels = "")
par(new=TRUE,fig=c(0,0.9,0.05,0.95))
nn=20
 
colors_inter <-rgb( ramp_inter(seq(0, 1, length = nn)), max = 255)
for(j in 1:nn){
	par(new=TRUE,fig=c(0,0.9,0.05,0.95))
	pie(y, col = colors_inter[j],radius = my_size*(inter_l+(inter_r-inter_l)*0.4*(1-j/nn)),xlab = "",ylab = "",border=0, labels = "")
}
par(new=TRUE,fig=c(0,0.9,0.05,0.95))
pie(1, col = shadow_color,radius = inter_l*my_size,xlab = "",ylab = "",border=0,labels = "")
par(new=TRUE,fig=c(0,0.9,0.05,0.95))
pie(1, col = "white",radius = (inter_l-(inter_r-inter_l)*0.15*narrow_tag)*my_size,xlab = "",ylab = "",border=0,labels = "")



par(new=TRUE,fig=c(0,0.9,0.05,0.95),font.main=4)
if (p==1){
	pie(y, col = ty,density = NULL, lty = NULL, angle = 45,main = paste(chrom_info[1], ":",ios_name_line, sep="")  ,radius = exon_r*my_size,xlab = "",ylab = "",border=0, labels = "",col.main="black",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0,)
}
else{
	pie(y, col = ty, angle = 45,radius = exon_r*my_size,xlab = "",ylab = "",border=0,labels = "",density = NULL,lty = NULL,clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)
}


par(new=TRUE,fig=c(0,0.9,0.05,0.95),font.main=4)
pie(y, col = ty_low,radius =(exon_r-(exon_r-exon_l)*0.1)*my_size,xlab = "",ylab = "",border=0, labels = "",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)

par(new=TRUE,fig=c(0,0.9,0.05,0.95))
pie(y, col = ty,radius =(exon_r-(exon_r-exon_l)*0.395)*my_size,xlab = "",ylab = "",border=0, labels = "",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)

nn=30
b3=ty

colors_exon<-rgb( ramp_exon(seq(0, 1, length = nn)), max = 255)
colors_intron<-rgb( ramp_intron(seq(0, 1, length = nn)), max = 255)
colors_UTR<-rgb( ramp_UTR(seq(0, 1, length = nn)), max = 255)
for(j in 1:nn){
	for(k in 1:length(ty)){
		if(ty_name[k] == "intron"){
			b3[k]=colors_intron[j]
			#cat(rainbow(n))
		}
		if(ty_name[k] == "CDS"){
			b3[k]=colors_exon[j]
			#cat(rainbow(n))
		}
		if(ty_name[k]=="five_prime_UTR" || ty_name[k]=="three_prime_UTR"){
			if(UTR_flag==1){
			   b3[k]=colors_UTR[j]
			 }
			 else{
                b3[k]=colors_exon[j]
			 }
			#cat(rainbow(n))
		}
	}
    par(new=TRUE,fig=c(0,0.9,0.05,0.95))
	pie(y, col = b3,radius = my_size*(exon_l+(exon_r-exon_l)*0.4*(1-j/nn)),xlab = "",ylab = "",border=0, labels = "",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)
}

b3=ty
for(k in 1:length(ty)){
		if(ty[k] != 0){
			b3[k]=shadow_color
			#cat(rainbow(n))
		}
}
par(new=TRUE,fig=c(0,0.9,0.05,0.95))
pie(y, col = b3,radius = exon_l*my_size,xlab = "",ylab = "",border=0,labels = "",clockwise=my_clockwise, init.angle = if(my_clockwise) 90 else 0)

par(new=TRUE,fig=c(0,0.9,0.05,0.95))
pie(x, col = "white",radius = (exon_l-(exon_r-exon_l)*0.05*narrow_tag)*my_size,xlab = "",ylab = "",border=0,labels = "")

}



par(new=TRUE,fig=c(0,1,0,1),font=1)
plot(1:100,1:100,type="n",xlab = "",ylab = "",axes = F)

text(d_w,d_h,chrom_info[2],pos=d_p1)
text(d_w,d_h,chrom_info[3],pos=d_p2)
x <- seq(x11,x12,length = 1000);
y <- seq(y11,y12,length = 1000);
lines(x, y)

mi_h=3;mi_m=1.2;sta_p=30;sta_l=80;m_l=4;m_shadow=0.2
if (exon_b==1){
	rect(sta_l+m_shadow, sta_p-m_shadow, sta_l+m_l+m_shadow, sta_p-mi_h-m_shadow,col =shadow_color,border=0,)
	rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =exon_color,border=0)
	if(UTR_flag==1){
		text(sta_l+m_l,sta_p-mi_h*0.6,"CDS", cex = 0.8,pos=4)
	}
	else{
		text(sta_l+m_l,sta_p-mi_h*0.6,"exon", cex = 0.8,pos=4)
	}
	sta_p=sta_p-mi_h-mi_m
}
if (intron_b==1){
	rect(sta_l+m_shadow, sta_p-m_shadow, sta_l+m_l+m_shadow, sta_p-mi_h-m_shadow,col =shadow_color,border=0,)
	rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =intron_color,border=0,)
	text(sta_l+m_l,sta_p-mi_h*0.6,"Intron", cex = 0.8,pos=4)
	sta_p=sta_p-mi_h-mi_m
}
if (inter_b==1){
	rect(sta_l+m_shadow, sta_p-m_shadow, sta_l+m_l+m_shadow, sta_p-mi_h-m_shadow,col =shadow_color,border=0,)
	rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =inter_color,border=0,)
	text(sta_l+m_l,sta_p-mi_h*0.6,"Intergenic", cex = 0.8,pos=4)
	sta_p=sta_p-mi_h-mi_m
}
if (UTR_b==1){
	rect(sta_l+m_shadow, sta_p-m_shadow, sta_l+m_l+m_shadow, sta_p-mi_h-m_shadow,col =shadow_color,border=0,)
	rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =UTR_color,border=0,)
	text(sta_l+m_l,sta_p-mi_h*0.6,"UTR", cex = 0.8,pos=4)
	sta_p=sta_p-mi_h-mi_m
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
gene_line_long=25
gene_line_wid=2
gene_line_posw=10
gene_line_posh=20
gene_line_interv=6


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
gene_line_posh<<-(gene_line_posh+gene_line_wid+gene_line_interv)
browers <- read.table(paste(file_dir1, "/",ios_name[p], sep="") ,stringsAsFactors = FALSE)

text(gene_line_posw-2,gene_line_posh+gene_line_wid+2,ios_name[p], cex = gene_line_cex,pos=4)

ty_name=browers[,1]
ty=browers[,1]
ty_low=browers[,1]

for(j in 1:length(ty)){
    if(ty[j]=="CDS"){
      ty[j]=exon_color
      ty_low[j]=exon_color_low
      exon_b=1
    }
#"#2E45F6"
    if(ty[j] == "intron"){
      ty[j]=intron_color
      ty_low[j]=intron_color_low
      intron_b=1
    }
    if(ty[j]=="five_prime_UTR" || ty[j]=="three_prime_UTR" || ty[j]=="UTR"){
      if(UTR_flag==1){
        ty[j]=UTR_color
          ty_low[j]=UTR_color_low
          UTR_b=1
      }
      else{
        ty[j]=exon_color
          ty_low[j]=exon_color_low
          exon_b=1
      }
    
    }
    if(ty[j] == "intergenic"){
      ty[j]="0"
      ty_low[j]="0"
      #cat(rainbow(n))
      inter_b=1
    }   
}


y=browers[,2]
x=1


colors_inter <-rgb( ramp_inter(seq(0, 1, length = nn)), max = 255)
for(j in 1:nn){
  #par(new=TRUE,fig=c(0,0.9,0.05,0.95))
  rect(gene_line_posw, gene_line_posh+gene_line_wid*0.15, gene_line_posw+gene_line_long,gene_line_posh+gene_line_wid*0.15+gene_line_wid*0.7*(1-j/nn),col =colors_inter[j],border=0,)
}



fun.geneline(y,ty,gene_line_posh+gene_line_wid*0.4,gene_line_posh+gene_line_wid)

nn=100
b3=ty

colors_exon<-rgb( ramp_exon(seq(0, 1, length = nn)), max = 255)
colors_intron<-rgb( ramp_intron(seq(0, 1, length = nn)), max = 255)
colors_UTR<-rgb( ramp_UTR(seq(0, 1, length = nn)), max = 255)
for(j in 1:nn){
  for(k in 1:length(ty)){
    if(ty_name[k] == "intron"){
      b3[k]=colors_intron[j]
      #cat(rainbow(n))
    }
    if(ty_name[k] == "CDS"){
      b3[k]=colors_exon[j]
      #cat(rainbow(n))
    }
    if(ty_name[k]=="five_prime_UTR" || ty_name[k]=="three_prime_UTR"){
      if(UTR_flag==1){
         b3[k]=colors_UTR[j]
       }
       else{
                b3[k]=colors_exon[j]
       }
      #cat(rainbow(n))
    }
  }
 fun.geneline(y,b3,gene_line_posh,gene_line_posh+gene_line_wid*0.5*(1-j/nn))
}

}

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

text(gene_line_posw,10,"gene struct", cex = 0.8,pos=4)
text(0,0,"(up to down,outer to inner)", cex = 0.7,pos=4)