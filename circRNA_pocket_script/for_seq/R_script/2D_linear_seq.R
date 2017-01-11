exon_b=0;intron_b=0;inter_b=0;UTR_b=0;show_gene_direction=FALSE;UTR_flag=1
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

my_size=0.8
seq_wid=0.05
seq_interv=0.05
show_mis=1
show_gap=1
mis_color=array("black",dim=c(2,1))
gap_color=array("white",dim=c(2,1))
seq_color1=array("grey",dim=c(2,1))
seq_color2=array("grey",dim=c(2,1))

w_vs_i=1;
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

fun.seqcolor <- function(col){
   if(col == "Blue") {seq_color1[1]<<-"blue";seq_color1[2]<<-"#0000A7"
        seq_color2[1]<<-"steelblue3";seq_color2[2]<<-"#3F76A4"
   }
    if(col == "Green") {seq_color1[1]<<-"green3";seq_color1[2]<<-"#00A400" 
        seq_color2[1]<<-"palegreen3";seq_color2[2]<<-"palegreen4"
   }
    if(col == "Yellow") {seq_color1[1]<<-"darkgoldenrod2";seq_color1[2]<<-"darkgoldenrod3"
       seq_color2[1]<<-"lightgoldenrod2";seq_color2[2]<<-"lightgoldenrod3"
   }
   if(col == "Red") {seq_color1[1] <<-"#EE1A1A";seq_color1[2] <<-"#CD0F0F"
       seq_color2[1] <<-"HotPink2";seq_color2[2] <<-"HotPink3"
   }  
   if(col == "Purple") {seq_color1[1]<<-"darkorchid3";seq_color1[2]<<-"darkorchid4"
       seq_color2[1]<<-"mediumpurple2";seq_color2[2]<<-"mediumpurple3"
   }    
   if(col == "Brown") {seq_color1[1]<<-"chocolate3";seq_color1[2]<<-"#A55217"
      seq_color2[1]<<-"#C1B26A";seq_color2[2]<<-"#A29453"
   }
    if(col == "Grey") { seq_color1[1]<<-"gray14";seq_color1[2]<<-"gray8"
       seq_color2[1]<<-"gray32";seq_color2[2]<<-"gray28";
   }
   if(col == "pink_blue") {seq_color1[1]<<-"rosybrown3";seq_color1[2]<<-"#AF8484"
       seq_color2[1]<<-"lightsteelblue3";seq_color2[2]<<-"#8A9AAF"
   }
    if(col == "red_gold") {seq_color1[1]<<-"brown2";seq_color1[2]<<-"brown3" 
     seq_color2[1]<<-"goldenrod2";seq_color2[2]<<-"goldenrod3"
   }
   if(col == "blue_green") {seq_color1[1]<<-"#1F1F99";seq_color1[2]<<-"darkblue"
     seq_color2[1]<<-"#1F761F";seq_color2[2]<<-"darkgreen"
   }
   if(col == "green_yellow") {seq_color1[1]<<-"#1F761F";seq_color1[2]<<-"darkgreen"
      seq_color2[1]<<-"yellow3";seq_color2[2]<<-"#B3B300"
   }
   if(col == "red_green") {seq_color1[1]<<-"#991F1F";seq_color1[2]<<-"darkred"
     seq_color2[1]<<-"#1F761F";seq_color2[2]<<-"darkgreen"
   }
    if(col == "blue_yellow") {seq_color1[1]<<-"#1F1F99";seq_color1[2]<<-"darkblue"
     seq_color2[1]<<-"yellow3";seq_color2[2]<<-"#B3B300"
   }
    if(col == "blue_red") {seq_color1[1]<<-"#1F1F99";seq_color1[2]<<-"darkblue"
     seq_color2[1]<<-"#991F1F";seq_color2[2]<<-"darkred"
   }
}

fun.width <- function(wid){
  line_cut<<-strsplit(wid,"_")[[1]]
  my_width=lapply(line_cut[1],as.numeric)[[1]]
  my_interval=lapply(line_cut[2],as.numeric)[[1]]
      w_vs_i<<-my_width/(my_interval+my_width)*2
      if(w_vs_i<0.4){
        w_vs_i<<-0.4
      }
       if(w_vs_i>1.2){
        w_vs_i<<-1.2
      }
}

Args <- commandArgs()
file_dir=Args[6]
file_format=Args[7]
file_output=Args[8]
n_line=1
con_format <- file(file_format, "r")
line=readLines(con_format,n=1)
show_gene_direction=FALSE

while(length(line) != 0 ) {

     if(n_line==1){
      fun.color(line)
     }
      if(n_line==2){
      fun.seqcolor(line)
     }
     if(n_line==3){
      circle_num_type<-strsplit(line,"\t")[[1]]
      fun.width(circle_num_type[2])
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


pdf(file_output)
#pdf("linear.pdf")
#file_dir_g="H:/d/circRNA_pocket_linux/output/Chr5_21715716_21724565/infog.txt"

#file_dir_g="H:/d/circRNA_pocket_linux/output/Chr5_3617186_3617467/infog.txt"
con_g <- file(paste(file_dir, "/infog.txt", sep=""), "r")
#con_g <- file(file_dir_g, "r")

par(fig=c(0,1,0,1),font=1)
plot(1:100,1:100,type="n",xlab = "",ylab = "",axes = F)
#plot(1:100,1:100,type="n",xlab = "",ylab = "")
#pdf("linear.pdf")
line=readLines(con_g,n=1)
pos_info<-strsplit(line,"\t")[[1]]
pos_chr=pos_info[1]
pos_sta=lapply(pos_info[2],as.numeric)[[1]]
pos_end=lapply(pos_info[3],as.numeric)[[1]]

line=readLines(con_g,n=1)
all_gene_range<-strsplit(line,"\t")[[1]]
all_star=lapply(all_gene_range[1],as.numeric)[[1]]
all_end=lapply(all_gene_range[2],as.numeric)[[1]]
all_len=all_end-all_star

p_sta<-(pos_sta-all_star)/all_len*100
p_end<-(pos_end-all_star)/all_len*100
rect(p_sta, 0, p_end, 100,col ="grey85",border=0,)
text(0,100,paste(pos_chr,":",pos_sta,"..",pos_end), cex = 0.7,pos=4);
line<-readLines(con_g,n=1) 
count_gene_trans<-strsplit(line,"\t")[[1]]
count_gene<-lapply(count_gene_trans[1],as.numeric)[[1]]
count_max_trans<-lapply(count_gene_trans[2],as.numeric)[[1]]
if(count_max_trans<6){
     count_max_trans<-6
}
for(i in 1:count_gene){
   
	line<-readLines(con_g,n=1)
	gene_info<-strsplit(line,"\t")[[1]]
	gene_name<-gene_info[1]
	count_trans<-lapply(gene_info[2],as.numeric)[[1]]
	gene_strd<-gene_info[3]
	g_hight<-50;
  g_range<-45;
  m_d_sta<-9999999999;
  m_d_end<-0;
	for(j in 1:count_trans){
		   line<-readLines(con_g,n=1)
           trans_info<-strsplit(line,"\t")[[1]]
           trans_name<-trans_info[1]
	       count_part<-lapply(trans_info[2],as.numeric)[[1]]
	       line<-readLines(con_g,n=1)
	       part_info<-strsplit(line,"\t")[[1]]
         
	       for(k in 1:count_part){
                draw_info<-strsplit(part_info[k],"-")[[1]]
                draw_type<-draw_info[1];
                draw_sta<-lapply(draw_info[2],as.numeric)[[1]]
                draw_end<-lapply(draw_info[3],as.numeric)[[1]]
                mcol<-UTR_color
                if(draw_type == "CDS"){
                   mcol<-exon_color;exon_b<<-1
                }
                if(draw_type == "intron"){
                   mcol<-intron_color;intron_b<<-1
                }
                if(draw_type == "UTR" || draw_type == "five_prime_UTR"|| draw_type == "three_prime_UTR"){
                     if(UTR_flag==1){
                         mcol<-UTR_color;UTR_b<<-1
                     }
                     else{
                         mcol<-exon_color;exon_b<<-1
                     }
                }
                d_sta<-(draw_sta-all_star)/all_len*100
                d_end<-(draw_end-all_star)/all_len*100
                m_w_vs_i=1
                #if(count_max_trans>1){
                m_w_vs_i<-w_vs_i
                #}
                #print(m_w_vs_i)
                rect(d_sta, g_hight, d_end, g_hight-g_range*0.5/count_max_trans*m_w_vs_i,col =mcol,border=0,)
                g_l<-g_hight-g_range*0.5/count_max_trans*m_w_vs_i*0.5;d_l<-1.5;
                g_inter=g_range*0.5/count_max_trans
                m_cex=0.6
                if(g_inter<3.5){
                    m_cex<-0.5
                }
                if(k==count_part && show_gene_direction && gene_strd == "+"){
                            x <- seq(d_end,d_end+d_l,length = 1000);
                            y <- seq(g_l,g_l,length = 1000);
                            lines(x, y)
                            x <- seq(d_end+0.65*d_l,d_end+d_l,length = 1000);
                            y <- seq(g_hight-g_range*0.1/count_max_trans*m_w_vs_i,g_l,length = 1000);
                            lines(x, y)
                }
                 if(k==1 && show_gene_direction && gene_strd == "-"){
                             x <- seq(d_sta,d_sta-d_l,length = 1000);
                             y <- seq(g_l,g_l,length = 1000);
                             lines(x, y)
                             x <- seq(d_sta-0.65*d_l,d_sta-d_l,length = 1000);
                             y <- seq(g_hight-g_range*0.1/count_max_trans*m_w_vs_i,g_l,length = 1000);
                             lines(x, y) 
                }
                if(i!=count_gene && k==1){
                             text(d_sta,g_hight+2,trans_name, cex = m_cex,pos=4);
                             
                }
                if(i==count_gene && k==count_part){
                             text(d_end,g_hight+2,trans_name, cex = m_cex,pos=2);
                             
                }
                if(d_sta<m_d_sta){
                                    m_d_sta<<-d_sta
                }
                if(d_end>m_d_end){
                                    m_d_end<<-d_end
                }

	       }
         g_hight<-(g_hight-g_range/count_max_trans)
	}
 
  if(i==count_gene){
       text(m_d_end,0,gene_name, cex = 0.7,pos=2,font=4)
  }
  else{
       if(i==1){
            text(m_d_sta,0,gene_name, cex = 0.7,pos=4,font=4)
        }else{
            text((m_d_sta+m_d_end)/2,0,gene_name, cex = 0.7,pos=3,font=4)
        }
       
  }
}

close(con_g)
m_axis=60
x <- seq(0,100,length = 1000);
y <- seq(m_axis,m_axis,length = 1000);
lines(x, y) 

fun.axis <- function(i,b){
      pos<-(i-all_star)/all_len*100
      x <- seq(pos,pos,length = 1000);
      y <- seq(m_axis,m_axis+1.5,length = 1000);
      lines(x, y) 
      t=i
      if(b==1){
          t<-paste(i/1000,"K")
      }
      text(pos,m_axis+0.5,t, cex = 0.6,pos=3)
}
for(i in all_star:all_end){
   if(all_len<=1000){
         if(i%%100 == 0){
                fun.axis(i,0)
         }
   }
   if(all_len>1000 && all_len<=2000){
         if(i%%200 == 0){
                fun.axis(i,0)
         }
   }
   if(all_len>2000 && all_len<=5000){
         if(i%%500 == 0){
                fun.axis(i,0)
         }
   }
    if(all_len>5000 && all_len<=10000){
         if(i%%1000 == 0){
                fun.axis(i,1)
         }
   }
   if(all_len>10000 && all_len<=20000){
         if(i%%2000 == 0){
                fun.axis(i,1)
         }
   }
   if(all_len>20000 && all_len<=50000){
         if(i%%5000 == 0){
                fun.axis(i,1)
         }
   }
   if(all_len>50000 && all_len<=100000){
         if(i%%10000 == 0){
                fun.axis(i,1)
         }
   }
   if(all_len>100000 && all_len<=200000){
         if(i%%20000 == 0){
                fun.axis(i,1)
         }
   }
   if(all_len>200000 && all_len<=500000){
         if(i%%50000 == 0){
                fun.axis(i,1)
         }
   }
    if(all_len>500000){
         if(i%%100000 == 0){
                fun.axis(i,1)
         }
   }
}



con_gene2 <- file(paste(file_dir, "/infog.txt", sep=""), "r")
line<-readLines(con_gene2,n=1)
site_arr<-strsplit(line,"\t")[[1]]
site_long<-(lapply(site_arr[3],as.numeric)[[1]]-lapply(site_arr[2],as.numeric)[[1]])
close(con_gene2)
seq_line_long=10
seq_line_wid=1
seq_line_posw=0
seq_line_posh=65
seq_line_interv=5
fun.seqdraw <- function(c,s,p,w){
     ss=array(0,dim=c(3,1))
    ss[1]=lapply(s[1],as.numeric)[[1]]
    ss[2]=lapply(s[2],as.numeric)[[1]]
    ss[3]=lapply(s[3],as.numeric)[[1]]

  
    mcol=array("",dim=c(3,0))
    mcol[2]=c[1]
   # pie(ss, col= mcol,radius = p,xlab = "",ylab = "",border=0,labels = "",clockwise=TRUE, init.angle = 90 )
    #rect((ss[1]+all_star)/all_len*100, seq_line_posh, (ss[1]+ss[2]+all_star)/all_len*100, seq_line_posh+seq_line_wid,col=c[1],border=0)
    rect((ss[1]+pos_sta-all_star)/all_len*100, seq_line_posh, (ss[1]+ss[2]+pos_sta-all_star)/all_len*100, seq_line_posh+seq_line_wid,col=c[1],border=0)
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
 #    x <- seq(pos+long,pos+long+2,length = 1000);
  #  y <- seq(seq_line_posh,seq_line_posh,length = 1000);
  #  lines(x, y)
    
         #x <- seq(pos+long+1,pos+long+2,length = 1000);
          #y <- seq(seq_line_posh+0.8,seq_line_posh,length = 1000);
          #lines(x, y)
   text(pos,seq_line_posh+seq_line_wid*0.5,s, cex =0.7,pos=4);
}

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

con_seq <- file(paste(file_dir, "/info_seql.txt", sep=""), "r")
seq_n<-lapply(readLines(con_seq,n=1) ,as.numeric)[[1]]
if(seq_n==1){
    seq_line_posh<-70
}
seq_line_posh<<-(seq_line_posh+(seq_line_wid+seq_line_interv)*(seq_n+1))
for(i in 1:seq_n){
   seq_line_posh<<-(seq_line_posh-(seq_line_wid+seq_line_interv))
   #seq_pos<<-(seq_pos+seq_wid+seq_interv)
   line<-readLines(con_seq,n=1) 
   text(seq_line_posw-2,seq_line_posh+seq_line_wid+2.5,line, cex = 0.6,pos=4)
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
    if(e_s>0){
          for(j in 1:e_s){
              fun.seq(seq_color2,2)     
           }
    }
   
}
close(con_seq)




con_seq <- file(paste(file_dir, "/info_seqc.txt", sep=""), "r")
#A=array(0,dim=c(1000,1))
#line<-readLines(con,n=1) 

seq_line_posh<-65
#plot(1:100,1:100,type="n",xlab = "",ylab = "")
#par(new=TRUE,fig=c(0,1,0,1),font=1)
seq_n<-lapply(readLines(con_seq,n=1) ,as.numeric)[[1]]
if(seq_n==1){
    seq_line_posh<-70
}
#seq_pos=my_size
#seq_pos<-(seq_pos+(seq_wid+seq_interv)*(seq_n+1))
seq_line_posh<<-(seq_line_posh+(seq_line_wid+seq_line_interv)*(seq_n+1))
for(i in 1:seq_n){
   seq_line_posh<<-(seq_line_posh-(seq_line_wid+seq_line_interv))
   #seq_pos<<-(seq_pos-(seq_wid+seq_interv))
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
    if(e_s>0){
          for(j in 1:e_s){
               fun.seq(seq_color2,1)     
          }
    }
  

}
close(con_seq)




mi_h=2;mi_m=1.2;sta_p=100;sta_l=90;m_l=2.8;m_shadow=0;m_cex=0.7
if (exon_b==1){
  rect(sta_l+m_shadow, sta_p-m_shadow, sta_l+m_l+m_shadow, sta_p-mi_h-m_shadow,col =shadow_color,border=0,)
  rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =exon_color,border=0)
  if(UTR_flag==1){
    text(sta_l+m_l,sta_p-mi_h*0.6,"CDS", cex = m_cex,pos=4)
  }
  else{
    text(sta_l+m_l,sta_p-mi_h*0.6,"exon", cex = m_cex,pos=4)
  }
  sta_p=sta_p-mi_h-mi_m
}
if (intron_b==1){
  rect(sta_l+m_shadow, sta_p-m_shadow, sta_l+m_l+m_shadow, sta_p-mi_h-m_shadow,col =shadow_color,border=0,)
  rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =intron_color,border=0,)
  text(sta_l+m_l,sta_p-mi_h*0.6,"Intron", cex = m_cex,pos=4)
  sta_p=sta_p-mi_h-mi_m
}
if (inter_b==1){
  rect(sta_l+m_shadow, sta_p-m_shadow, sta_l+m_l+m_shadow, sta_p-mi_h-m_shadow,col =shadow_color,border=0,)
  rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =inter_color,border=0,)
  text(sta_l+m_l,sta_p-mi_h*0.6,"Intergenic", cex = m_cex,pos=4)
  sta_p=sta_p-mi_h-mi_m
}
if (UTR_b==1){
  rect(sta_l+m_shadow, sta_p-m_shadow, sta_l+m_l+m_shadow, sta_p-mi_h-m_shadow,col =shadow_color,border=0,)
  rect(sta_l, sta_p, sta_l+m_l, sta_p-mi_h,col =UTR_color,border=0,)
  text(sta_l+m_l,sta_p-mi_h*0.6,"UTR", cex = m_cex,pos=4)
  sta_p=sta_p-mi_h-mi_m
}
if(show_gene_direction ){
    x <- seq(sta_l+m_l,sta_l,length = 1000);
    y <- seq(sta_p-mi_h,sta_p-mi_h,length = 1000);
    lines(x, y)

         x <- seq(sta_l+0.7*m_l,sta_l+m_l,length = 1000);
         y <- seq(sta_p+0.5*mi_h-mi_h,sta_p-mi_h,length = 1000);
         lines(x, y)
    text(sta_l+m_l,sta_p-mi_h*0.6,"Direction", cex = m_cex,pos=4)
}