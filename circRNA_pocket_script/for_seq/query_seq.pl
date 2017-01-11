#!/usr/bin/perl -w
use strict;
use File::Basename;
my $my_dir=dirname(__FILE__);
my $input_config =$my_dir."/config.txt";

#my $input1 = "seq_blast.out";
#my $input2 = "seq.fa";
#my $input_anno = "H:/d/annote/all.gff3";
#my $input_anno = "F:/d/annote/rice_refFlat_MSU.txt";
#my $anno_type="gff3";
#my $input_refer = "H:/d/soft_test_small/rice_data/all.fa";
#open IN1, '<'.$input1;
#open IN2, '<'.$input2;
#open IN_anno, '<'.$input_anno;
#open IN_refer, '<'.$input_refer;

#open IN_config, '<'.$input_config or die "Can't open $input_config. $!";
#opendir(IN_config,$input_config) or die "Can't open $input_config. $!";
my %h;
my @qseqid;my @sseqid;my @pident;
my @leng;my @mismatch;my @gapopen;
my @qstart;my @qend;my @sstart;
my @send;my @evalue;my @bitscore;
my @qseq;my @sseq;my @strd;
my $count=0;
my $maxdis=5000;
my %h_seq;my $s4;
my %h_anno_up;
my %h_anno_down;


my %chrom_gene;
my %gene_mRNA;
my %gene_up;
my %gene_down;
my %gene_strd;
my %m_trans;

my $anno_type;
my $input_qs;
my $input_blast;
my $input_anno;
my $input_blastdb;
my $color_type="3D";
my $output_dir="circRNA_pocket_seq_output";

my $b_test_config=test_config();
my $b_test_format=test_format();
if($b_test_config!=1){die "Exit.";}
if($b_test_format!=1){die "Exit.";}
sub test_config
{
    open IN_config, '<'.$input_config or die "Can't open $input_config. $!";
    my @c_name=("blast_database","annotation_file","annotation_type","input_query_sequence_file","color_type","output_dir","max_dis");
    my @b_name=(0,0,0,0,0,0,0);
    while (<IN_config>){
          my $s=$_;
          chomp($s);
          my @c=split(/=/,$s);
          if (substr($s,0,1) eq "#" || $s eq "") {next;}
     
          if($c[0] eq "blast_database"){
             $b_name[0]=1;
             $input_blastdb=$c[1];next;
          }
          if($c[0] eq "annotation_file"){
             $b_name[1]=1;
             $input_anno=$c[1];next;
          }
          if($c[0] eq "annotation_type"){
             $b_name[2]=1;
             $anno_type=$c[1];next;
          }
          if($c[0] eq "input_query_sequence_file"){
             $b_name[3]=1;
             $input_qs = $c[1];next;
          }
          if($c[0] eq "color_type"){
             $b_name[4]=1;
             $color_type = $c[1];next;
          }
          if($c[0] eq "output_dir"){
             $b_name[5]=1;
             $output_dir = $c[1];next;
          }
          if($c[0] eq "max_dis"){
             $b_name[6]=1;
             $maxdis = $c[1];next;
          }
    }
    close IN_config;
    for my $i(0..3){
         if($b_name[$i]==0){
            print "$input_config. Can not find the item $c_name[$i].\n";
            return 0;
         }
    }
    if($maxdis<100){
         print "$input_config. max_dis should be at least 100.\n";
         return 0;
    }
    return 1;
}
sub test_format
{
     #cannot open file 'output/Chr1_39702072_39702154/\info_seqc.txt': No such file or directory
     my $input_format =$my_dir."/format.txt";
     open IN_format, '<'.$input_format or die "Can't open $input_format. $!";
     my $count=0;
      my $x1;my $y1;my $x2;my $y2;my $x3;my $y3;my $x4;
     while (<IN_format>){
          my $s=$_;
           chomp($s);
           $count++;
           if($count==1){
               if($s eq "Blue" || $s eq "Red" || $s eq "Yellow" || $s eq "Green" || $s eq "Purple" || $s eq "Brown" || $s eq "Grey" || $s eq "pink_white_blue" || $s eq "red_wood_gold" || $s eq "pink_white_blue" || $s eq "red_wood_gold" || $s eq "blue_red_green" || $s eq "blue_yellow_green")
               {;}
               else{
                  print "$input_format\tline 1.The color of gene should be (Blue;Red;Yellow;Green;Purple;Brown;Grey;pink_white_blue;red_wood_gold;pink_white_blue;red_wood_gold;blue_red_green;blue_yellow_green)\n";
                  return 0;
               }
               next;
           }
           if($count==2){
               if($s eq "Blue" || $s eq "Green" || $s eq "Yellow" || $s eq "Red" || $s eq "Purple" || $s eq "Brown" || $s eq "Grey" || $s eq "pink_blue" || $s eq "red_gold" || $s eq "blue_green" || $s eq "green_yellow" || $s eq "red_green" || $s eq "blue_yellow")
               {;}
               else{
                   print "$input_format\tline 2.The color of input sequences should be (Blue;Green;Yellow;Red;Purple;Brown;Grey;pink_blue;red_gold;blue_green;green_yellow;red_green;blue_yellow)\n";
                   return 0;
               }
               next;
           }
          if($count==3){
              my @c=split(/\t/,$s);
              my @c1=split(/_/,$c[0]);my @c2=split(/_/,$c[1]);
              $x1=$c1[0];$y1=$c1[1];$x2=$c2[0];$y2=$c2[1];
              if($x1<2 || $x2<2 || $x1>10 || $x2>10 || $y2<0 || $y2>10 || $y1<2 || $y1>6){
                  print "$input_format\tline 3.x1_y1 x2_y2 should be (3 <= x1,x2 <= 10; 0 <= y2 <= 10;2 <= y1 <=6)\n";
                  return 0;
              }
              next;
          }
          if($count==4){
              my @c=split(/_/,$s);
              $x3=$c[0];$y3=$c[1];
              if($x3<1 || $x3>5 || $y1<3 || $y1>5){
                  print "$input_format\tline 4.x3_y3 should be (1 <= x3 <= 5; 1 <= y3 <= 5)\n";
                  return 0;
              }
              next;
          }
          if($count==5){
              if((not $s eq "show_UTR") && (not $s eq "no_UTR")){
                  print "$input_format\tline 5.It should be (show_UTR;no_UTR)\n";return 0;
              }
              next;
          }
           if($count==6){
              if((not $s eq "show_gene_direction") && (not $s eq "no_gene_direction")){
                  print "$input_format\tline 6.It should be (show_gene_direction;no_gene_direction)\n";return 0;
              }
              next;
          }
          if($count==7){
              if((not $s eq "show_mis") && (not $s eq "no_mis")){
                  print "$input_format\tline 7.It should be (show_mis;no_mis)\n";return 0;
              }
              next;
          }
          if($count==8){
              if((not $s eq "show_gap") && (not $s eq "no_gap")){
                  print "$input_format\tline 8.It should be (show_gap;no_gap)\n";return 0;
              }
              next;
          }
          if($count==9){
              $x4=$s;
              if($s<0.5 && $s> 0.9){
                  print "$input_format\tline 9.x4,The radius of gene in the circle format should be (0.5 <= x4 <= 0.9)\n";return 0;
              }
              next;
          }
     }
     close IN_format;
     return 1;
}


mkdir($output_dir);
my $command="";
$command=qq(blastn -db  $input_blastdb -query  $input_qs -out $output_dir/seq_blast.out -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" -word_size 20 -penalty -5 -gapopen 15 -gapextend 3);
my $b_command_res=system($command);
if($b_command_res!=0){
   die "blastn failed!";
}
$input_blast=$output_dir."/seq_blast.out";

open IN_qs, '<'.$input_qs or die "Can't open $input_qs. $!";
open IN_blast, '<'.$input_blast or die "Can't open $input_blast. $!";
open IN_anno, '<'.$input_anno or die "Can't open $input_anno. $!";
open IN_refer, '<'.$input_blastdb or die "Can't open $input_blastdb. $!";

print "Reading genome annotation file..\n";
if ($anno_type eq "gff3") {get_anno_gff3();}
if ($anno_type eq "gtf") {get_anno_gtf();}
if ($anno_type eq "refFlat") {get_anno_refFlat();}


sub get_anno_gff3
{

 while (<IN_anno>){
      my $s=$_;
     chomp($s);
     if (substr($s,0,1) eq "#") {next;}
     if ($s eq "") {next;}
     
    my @c=split(/\t/,$s);
    my $chrom=$c[0];
    my $type=$c[2];
    my $sta=$c[3]-1;my $end=$c[4];
    if ($type eq "exon"){
        $h_anno_up{$c[0]}{$sta}=1;
        $h_anno_down{$c[0]}{$end}=1;
    }
    if ($type eq "gene"){
        my $ID=get_it_gff3("ID",$c[8]);
        push(@{$chrom_gene{$chrom}},$ID);
        $gene_up{$ID}=$sta;
        $gene_down{$ID}=$end;
        $gene_strd{$ID}=$c[6];
    }
    if (($type eq "five_prime_UTR")||($type eq "CDS")||($type eq "three_prime_UTR")) {
        my $Parent=get_it_gff3("Parent",$c[8]);
        push(@{$m_trans{$Parent}},$type."-".$sta."-".$end);
    }
    if ($type eq "mRNA") {
        my $Parent=get_it_gff3("Parent",$c[8]);
        my $ID=get_it_gff3("ID",$c[8]);
        push(@{$gene_mRNA{$Parent}},$ID);
    }
    
  }
 close IN_anno;
 sub get_it_gff3
 {
     my($gt,$s)=@_;
     my @c=split(/;/,$s);
     for my $i(@c){
          my($l,$r)=split(/=/,$i);
          if ($l eq $gt) {return $r; }       
     }
     return;
  }
=for
  my $kkk="Chr5_21715716_21724565";
  mkdir("output/".$kkk);
  print_gene_struct($kkk);
  my $aaa=1;

   my $kkk="Chr5_21715716_21724565";
   mkdir("output/".$kkk);
   print_gene_beside($kkk);
   my $aaa=1;
=cut
}

for my $k(keys %gene_mRNA){
    my $b=0;
    if(not exists($gene_up{$k})){next;}
    for my $i(@{$gene_mRNA{$k}}){
        if(exists($m_trans{$i})){$b=1;next;}
        if($b==1){last;}
    }
    

    if($b==0){
         for my $i(@{$gene_mRNA{$k}}){
            my $gsta=$gene_up{$k};
            my $gend=$gene_down{$k};
           # print "myprint k:".$k."\n";
            push(@{$m_trans{$i}},"UTR-".$gsta."-".$gend);
       }
    }
}

sub get_anno_gtf
{

   my %gene_mRNA_k;my %m_trans_temp;
   while (<IN_anno>){
     my $s=$_;
     chomp($s);
     if (substr($s,0,1) eq "#") {next;}
     if ($s eq "") {next;}
     
     my @c=split(/\t/,$s);
     my $chrom=$c[0];
     my $type=$c[2];
     my $sta=$c[3]-1;my $end=$c[4];
     if ($type eq "exon") {
         $h_anno_up{$c[0]}{$sta}=1;
         $h_anno_down{$c[0]}{$end}=1;
     }
    
     if (($type eq "exon")||($type eq "CDS")) {
        my $transcript_id=get_it_gtf("transcript_id",$c[8]);
        push(@{$m_trans_temp{$transcript_id}},$type."-".$sta."-".$end);
        my $gene_id=get_it_gtf("gene_id",$c[8]);
        $gene_mRNA_k{$chrom.";".$gene_id.";".$c[6]}{$transcript_id}=1;
     }
  }
  close IN_anno;
  for my $k(keys %m_trans_temp)
  {
     my $max_cds=0;my $min_cds=9999999999;
     my %h_cds_up; my %h_cds_down;
     for my $i(@{$m_trans_temp{$k}}){
          my ($type,$sta,$end)=split(/-/,$i);
          if ($type eq "CDS") {
               $h_cds_up{$sta}=1;$h_cds_down{$end}=1;
              if ($max_cds<$end) {$max_cds=$end;}
              if ($min_cds>$sta) {$min_cds=$sta;}
          }
     }

     for my $i(@{$m_trans_temp{$k}}){
          my ($type,$sta,$end)=split(/-/,$i);
          if ($type eq "exon") {
               if (exists($h_cds_up{$sta}) && exists($h_cds_down{$end})) {
                    push(@{$m_trans{$k}},"CDS-".$sta."-".$end);next;
               }
               if (exists($h_cds_up{$sta}) &&(not exists($h_cds_down{$end}))) {
                    push(@{$m_trans{$k}},"CDS-".$sta."-".$max_cds);
                    push(@{$m_trans{$k}},"UTR-".$max_cds."-".$end);next;
               }
               if ((not exists($h_cds_up{$sta})) &&( exists($h_cds_down{$end}))) {
                    push(@{$m_trans{$k}},"UTR-".$sta."-".$min_cds);
                    push(@{$m_trans{$k}},"CDS-".$min_cds."-".$end);next;
               }
               if ((not exists($h_cds_up{$sta})) &&(not exists($h_cds_down{$end}))) {
                    if ($min_cds>$end || $max_cds< $sta) {push(@{$m_trans{$k}},"UTR-".$sta."-".$end);next; }
                    if ($min_cds>$sta || $max_cds< $end) {
                        push(@{$m_trans{$k}},"UTR-".$sta."-".$min_cds);
                        push(@{$m_trans{$k}},"CDS-".$min_cds."-".$max_cds);
                        push(@{$m_trans{$k}},"UTR-".$max_cds."-".$end);
                    }
               }
          }
     }
  }
  for my $k_g(keys %gene_mRNA_k){
      my ($chr,$g_id,$std)=split(/;/,$k_g);
      push(@{$chrom_gene{$chr}},$g_id);
      $gene_strd{$g_id}=$std;
      my $end=0;my $sta=9999999999;
      for my $k_t(keys %{$gene_mRNA_k{$k_g}}){
           push(@{$gene_mRNA{$g_id}},$k_t);
           for my $i(@{$m_trans{$k_t}}){
                my($type,$msta,$mend)=split(/-/,$i);
                if ($msta<$sta) {$sta=$msta;}
                if ($mend>$end) {$end=$mend;}
           }
      }
      $gene_up{$g_id}=$sta;
      $gene_down{$g_id}=$end; 
  }
  sub get_it_gtf
  {
     my($gt,$s)=@_;
     my @c=split(/ /,$s);
     for (my $i=0;$i<$#c;$i+=2)
     {
        if ($c[$i] eq $gt) {
          return substr($c[$i+1],1,length($c[$i+1])-3);
        }
     }
     return;
  }
 
}



sub get_anno_refFlat
{
   my %chrom_gene_temp;
   while (<IN_anno>){
      my $s=$_;
      chomp($s);
      if (substr($s,0,1) eq "#") {next;}
      if ($s eq "") {next;}
      my @c=split(/\t/,$s);
      my $tsta=$c[4];my $tend=$c[5];
      my $gene=$c[0];my $trans=$c[1];
      my $strd=$c[3];
      my $cds_l=$c[6]; my $cds_r=$c[7];
      
      if (not exists($gene_up{$gene})) {$gene_up{$gene}=$tsta;}
      else{
          if ($tsta<$gene_up{$gene}) {$gene_up{$gene}=$tsta;}
      }
      if (not exists($gene_down{$gene})) {$gene_down{$gene}=$tend;}
      else{
          if ($tend>$gene_down{$gene}) {$gene_down{$gene}=$tend;}
      }
     
      
      $chrom_gene_temp{$c[2].";".$gene}=1;
      $gene_strd{$gene}=$strd;
      push(@{$gene_mRNA{$gene}},$trans);
      my @c_l=split(/,/,$c[9]);
      my @c_r=split(/,/,$c[10]);
      for my $i(0..$#c_l){
          my $sta=$c_l[$i];my $end=$c_r[$i];
          if ($sta eq "" and $end eq ""){next;}
          $h_anno_up{$gene}{$sta}=1;
          $h_anno_down{$gene}{$end}=1;
          if ($sta>= $cds_l && $end <= $cds_r) {push(@{$m_trans{$trans}},"CDS-".$sta."-".$end);next;}
          if ($sta>= $cds_r || $end <= $cds_l) {push(@{$m_trans{$trans}},"UTR-".$sta."-".$end);next;}
          if ($sta<= $cds_r && $sta >= $cds_l && $end > $cds_r) {
               push(@{$m_trans{$trans}},"CDS-".$sta."-".$cds_r);
               push(@{$m_trans{$trans}},"UTR-".$cds_r."-".$end);next;
          }
          if ($end<= $cds_r && $end >= $cds_l && $sta < $cds_l) {
               push(@{$m_trans{$trans}},"UTR-".$sta."-".$cds_l);
               push(@{$m_trans{$trans}},"CDS-".$cds_l."-".$end);next;
          }
          if ($sta < $cds_l && $end > $cds_r) {
               push(@{$m_trans{$trans}},"UTR-".$sta."-".$cds_l);
               push(@{$m_trans{$trans}},"CDS-".$cds_l."-".$cds_r);
               push(@{$m_trans{$trans}},"UTR-".$cds_r."-".$end);next;
          }
      }
    }
    close IN_anno;
   for my $k(keys %chrom_gene_temp){
        my @c_g=split(/;/,$k);
        push(@{$chrom_gene{$c_g[0]}},$c_g[1]);
   }
 
}

print "Reading genome file:\n";
my %h_refer;
while (<IN_refer>){
      my $s=$_;
     chomp($s);
     my $s3=substr($s,0,1);
    if ($s3 eq ">") {
       my @c=split(/ /,$s);
       $s4=substr($c[0],1);
       print $s4."\n";
    }
    else
    {    
       $h_refer{$s4}.=$s;      
    }   
}
close IN_refer;

print "Reading query sequence file:\n";
while (<IN_qs>){
     my $s=$_;
     chomp($s);
     my $s3=substr($s,0,1);
     if ($s3 eq ">") {
        my @c=split(/ /,$s);
        $s4=substr($c[0],1);
        print $s4."\n";
     }
     else{    
        $h_seq{$s4}.=$s;      
     }   
}
close IN_qs;

my %h_seqs;
for my $k (keys %h_seq){
     $h_seqs{$k}=$h_seq{$k};
     $h_seq{$k}=length($h_seq{$k});
}

while (<IN_blast>){
     my $s=$_;
     chomp($s);
     my @c=split(/,/,$s);
     if ($c[10]>1.00e-005) {next;}
     $c[6]--;
     if ($c[9]>$c[8]) {push(@strd,1);}
     else{
          
          my $seq_l = $h_seq{$c[0]};
	      my $temp ;
          $temp = ($seq_l - $c[6]) ;
	      $c[6] = ($seq_l - $c[7]) ;
	      $c[7] = $temp;
	 # $temp = $c[8];
	 # $c[8] = $c[9];
	#  $c[9] = $temp;
          ($c[8],$c[9])=($c[9],$c[8]);
          push(@strd,0);
     }
    $c[8]--;
    push(@qseqid,$c[0]);
    push(@sseqid,$c[1]);
    push(@pident,$c[2]);
    push(@leng,$c[3]);
    push(@mismatch,$c[4]);
    push(@gapopen,$c[5]);
    push(@qstart,$c[6]);
    push(@qend,$c[7]);
    push(@sstart,$c[8]);
    push(@send,$c[9]);
    push(@evalue,$c[10]);
    push(@bitscore,$c[11]);
    push(@qseq,$c[12]);
    push(@sseq,$c[13]);
    
    
    $count++;
}
close IN_blast;

my @order;
for my $i(0..($count-1)){
     push(@order,$i);
}
my %spl;
my %spr;
my $end_ad=7;
for my $i(sort{$qstart[$a]<=>$qstart[$b]} @order){
     if ($qstart[$i]<$end_ad) {
         if ($qstart[$i]==0 || ($evalue[$i]<1.00e-007)) {
            push(@{$spl{$qseqid[$i].";".$sseqid[$i].";".$strd[$i].";".$qstart[$i].";".$qend[$i].";".$sstart[$i].";".$send[$i]}},$i);
          }
     }
     else{
          for my $k(keys %spl){
               my($qid,$sid,$srd,$qsta,$qend,$ssta,$sed)=split(/;/,$k);
               if ((not $qid eq $qseqid[$i])|| $srd!=$strd[$i] || (not $sid eq $sseqid[$i])) {next;}            
               if ($sed < $sstart[$i] && $qend - $qstart[$i]<10 && abs($sed - $sstart[$i]) < $maxdis) {
                   push(@{$spl{$k}},$i);
               }              
          }
     }
}

for my $i(sort{$qend[$b]<=>$qend[$a]} @order){
     my $seql=$h_seq{$qseqid[$i]};
     if(abs($qend[$i] - $seql) < $end_ad){
          if (abs($qend[$i] - $seql)==0 || ($evalue[$i]<1.00e-007)) {
            push(@{$spr{$qseqid[$i].";".$sseqid[$i].";".$strd[$i].";".$qstart[$i].";".$qend[$i].";".$sstart[$i].";".$send[$i]}},$i);
          }
     }
     else{
          for my $k(keys %spr){
               my($qid,$sid,$srd,$qsta,$qend,$ssta,$sed)=split(/;/,$k);
               if ((not $qid eq $qseqid[$i]) || $srd!=$strd[$i] || (not $sid eq $sseqid[$i])) {next;}            
               if ($ssta > $send[$i] && $qend[$i] - $qsta<10 &&abs($send[$i] - $ssta) < $maxdis) {
                   push(@{$spr{$k}},$i);
               }              
          }
     }
}
my $spl_b=0;
while ($spl_b==0) {
     my $b_out=0;
     for my $k(keys %spl){
         for my $i(1..$#{$spl{$k}}){

              my $n1=${$spl{$k}}[$i-1];
              my $n2=${$spl{$k}}[$i];
              if (($qend[$n1]-$qstart[$n2] >10) || ($sstart[$n2] < $send[$n1])){
                 
                   my @a1=@{$spl{$k}};my @a2=@{$spl{$k}};
                   delete($spl{$k});
                   @{$spl{$k."_1"}}=splice(@a1,$i-1,1);
                   @{$spl{$k."_2"}}=splice(@a1,$i,1);
                   $b_out=1;
                   next;
              }
          }
          if ($b_out==1) {next;}        
     }
     if ($b_out==0) {$spl_b=1;}
}

my $spr_b=0;
while ($spr_b==0) {
     my $b_out=0;
     for my $k(keys %spr){
         for my $i(1..$#{$spr{$k}}){

              my $n1=${$spr{$k}}[$i-1];
              my $n2=${$spr{$k}}[$i];
              if (($qend[$n2]-$qstart[$n1] >10) || ($sstart[$n1] < $send[$n2])){
                 
                   my @a1=@{$spr{$k}};my @a2=@{$spr{$k}};
                   delete($spr{$k});
                   @{$spr{$k."_1"}}=splice(@a1,$i-1,1);
                   @{$spr{$k."_2"}}=splice(@a1,$i,1);
                   $b_out=1;
                   next;
              }
          }
         if ($b_out==1) {next;}        
     }
     if ($b_out==0) {$spr_b=1;}
}

my @drawl;
my @drawr;
my @drawl_ad;
my @drawr_ad;
my %draw_pos;
my %draw_linear_pos;
#my @drawlinear;
my @drawlinear;
my @drawlinear_ad;
my $count_draw=0;
my $count_linear_draw=0;
my %query_type;
for my $kl(keys %spl){
     my @cl=split(/;/,$kl);
      my $nls=${$spl{$kl}}[0];
      #my $nle=${$spl{$kl}}[-1];
      #my $suml=getsum(@{$spl{$kl}});
     for my $kr(keys %spr){
          my @cr=split(/;/,$kr);
          #my $sumr;my $suml;
          if (($cl[0] eq $cr[0]) && ($cl[1] eq $cr[1]) && $cl[2]==$cr[2]) {
               my $nrs=${$spr{$kr}}[0];
               #my $nre=${$spr{$kr}}[-1];
              # if($send[$nrs]< $sstart[$nls] || abs($send[$nrs]-$sstart[$nls])>$maxdis){next; }
              if($send[$nrs]> $sstart[$nls] || abs($send[$nrs]-$sstart[$nls])>$maxdis){next; }
              #my @spa=@{$spl{$kl}};my @spb=@{$spr{$kr}};
              match_it(\@{$spl{$kl}},\@{$spr{$kr}});
=for
               $sumr=getsum(@{$spr{$kr}});
               if ($suml+$sumr>=0.9) {
                    push(@drawl,[@{$spl{$kl}}]);
                    push(@drawr,[@{$spr{$kr}}]);
                    my @a1=get_bound_ad_l(@{$spl{$kl}});
                    my @a2=get_bound_ad_r(@{$spr{$kr}});
                    my ($a11,$a21)=check_splice_site(\@a1,\@a2,$cl[1]);
                    push(@drawl_ad,[@$a11]);
                    push(@drawr_ad,[@$a21]);
               }
=cut
          }
     }
}

for my $kl(keys %spl){
      
      my @c=@{$spl{$kl}};
      #my @c1=get_bound_ad_l(\@c,0);
      my $suml=getsum(@c);
      if ($suml>0.9) {
        my $bits=0;
        for my $i(@c){$bits+=$bitscore[$i];}
       # $bits/=($#c+1);
        my $nls=${$spl{$kl}}[0];
         my $nle=${$spl{$kl}}[-1];
         my $qid=$qseqid[$nls];
         my $pos=$sseqid[$nls].":".$sstart[$nls]."..".$send[$nle];
         my $k=$sseqid[$nls]."_".$sstart[$nls]."_".$send[$nle];
         push(@{$query_type{$qid}},"linear_pos,".$pos.",".$bits);
         push(@drawlinear,[@{$spl{$kl}}]);
         my @a1=get_bound_ad_l(\@c,0);
          push(@drawlinear_ad, [@a1]);
          push(@{$draw_linear_pos{$k}},$count_linear_draw);
            $count_linear_draw++;
      }
}

for my $k(keys %draw_linear_pos){
      mkdir($output_dir."/".$k);
    # print_gene_struct($k);
     print_gene_beside($k);
     open OUT, '> '.$output_dir.'/'.$k."/info_seqc.txt";
     my $s=($#{$draw_linear_pos{$k}}+1)."\n";
     my ($chr,$sta,$end)=split(/_/,$k);
     for my $i(@{$draw_linear_pos{$k}}){
          my @al=@{$drawlinear[$i]};
          my $std="+";
          if ($strd[$al[0]]==0) {$std="-";} 
          $s.=$qseqid[$al[0]]."\n".$std."\n".($#al+1)."\t0\n";
          my $s1;
          my @l=@{$drawlinear_ad[$i]};
          for my $j(0..$#l){
               my($tem1,$tem2,$jl,$jr)=split(/_/,$l[$j]);
               $s1.=($jl-$sta)."\t".($jr-$jl)."\t".($end-$jr)."\n".get_misgap_s($al[$j],$sta,$end)."\n";
          }
          $s.=$s1;
     }
     print OUT $s;
     close OUT;
     open OUT, '> '.$output_dir.'/'.$k."/info_seql.txt";
     $s=($#{$draw_linear_pos{$k}}+1)."\n";
     for my $i(@{$draw_linear_pos{$k}}){
          my @al=@{$drawlinear[$i]};
          my $qid=$qseqid[$al[0]];
          my $std="+";
          if ($strd[$al[0]]==0) {$std="-";} 
          $s.=$qid."\n".$h_seq{$qid}."\n".$std."\n".($#al+1)."\t0\n";
          my $s1;
          my @l=@{$drawlinear_ad[$i]};
          for my $j(0..$#l){
               my($jl,$jr,$tem1,$tem2)=split(/_/,$l[$j]);
               $s1.=$jl."\t".($jr-$jl)."\n".get_misgap_p($al[$j])."\n";
          }
          $s.=$s1;
     }
     print OUT $s;
     close OUT;
     
    draw_linear_R($output_dir."/".$k.'/');
}
sub match_it
{
     my ($sa,$sb)=@_;
      my @spa=@$sa;my @spb=@$sb;
       my @pa=@$sa;my @pb=@$sb;
       if ($#spa<0 || $#spb<0) {return;}
     
     my $b=0;
      for my $i(0..$#spa){
          if ($i>0) { pop(@pa);}
          @pb=@$sb;
          for my $j(0..$#spb){
               if ($j>0) {pop(@pb);}
               if ($qend[$pa[-1]] <$qstart[$pb[-1]]) {last;}
               if (abs($qend[$pa[-1]] - $qstart[$pb[-1]])<10) {
                    match_it_ad(\@pa,\@pb);
                    $b=1;last;
               }
          }
          if ($b==1) {last;}  
      }
}


sub match_it_ad
{
     my ($sa_ad,$sb_ad)=@_;
     my @c1=@$sa_ad;my @c2=@$sb_ad;
     push(@drawl,[@c1]);
     push(@drawr,[@c2]);
     my $minu=$qend[$c1[-1]] - $qstart[$c2[-1]];
     my @a1=get_bound_ad_l(\@c1,$minu);
     my @a2=get_bound_ad_r(@c2);
     my ($a11,$a21)=check_splice_site(\@a1,\@a2,$sseqid[$c1[0]]);
     push(@drawl_ad,[@$a11]);
     push(@drawr_ad,[@$a21]);
}



for my $k(keys %draw_pos){
     mkdir($output_dir."/".$k);
    # print_gene_struct($k);
     print_gene_struct($k);
     print_gene_beside($k);
     open OUT, '>'.$output_dir.'/'.$k."/info_seqc.txt";
     my $s=($#{$draw_pos{$k}}+1)."\n";
     my ($chr,$sta,$end)=split(/_/,$k);
     for my $i(@{$draw_pos{$k}}){
          my @al=@{$drawl[$i]};my @ar=@{$drawr[$i]};
          my $std="+";
          if ($strd[$al[0]]==0) {$std="-";} 
          $s.=$qseqid[$al[0]]."\n".$std."\n".($#al+1)."\t".($#ar+1)."\n";
          my $s1;
          my @l=@{$drawl_ad[$i]};my @r=@{$drawr_ad[$i]};
          for my $j(0..$#l){
               my($tem1,$tem2,$jl,$jr)=split(/_/,$l[$j]);
               $s1.=($jl-$sta)."\t".($jr-$jl)."\t".($end-$jr)."\n".get_misgap_s($al[$j],$sta,$end)."\n";
          }
          for my $j(0..$#r){
               my($tem1,$tem2,$jl,$jr)=split(/_/,$r[$j]);
               $s1.=($jl-$sta)."\t".($jr-$jl)."\t".($end-$jr)."\n".get_misgap_s($ar[$j],$sta,$end)."\n";
          }
          $s.=$s1;
     }
     print OUT $s;
     close OUT;
     open OUT, '> '.$output_dir.'/'.$k."/info_seql.txt";
     $s=($#{$draw_pos{$k}}+1)."\n";
     for my $i(@{$draw_pos{$k}}){
          my @al=@{$drawl[$i]};my @ar=@{$drawr[$i]};
          my $qid=$qseqid[$al[0]];
          my $std="+";
          if ($strd[$al[0]]==0) {$std="-";} 
          $s.=$qid."\n".$h_seq{$qid}."\n".$std."\n".($#al+1)."\t".($#ar+1)."\n";
          my $s1;
          my @l=@{$drawl_ad[$i]};my @r=@{$drawr_ad[$i]};
          for my $j(0..$#l){
               my($jl,$jr,$tem1,$tem2)=split(/_/,$l[$j]);
               $s1.=$jl."\t".($jr-$jl)."\n".get_misgap_p($al[$j])."\n";
          }
          for my $j(0..$#r){
               my($jl,$jr,$tem1,$tem2)=split(/_/,$r[$j]);
               $s1.=$jl."\t".($jr-$jl)."\n".get_misgap_p($ar[$j])."\n";  
          }
          $s.=$s1;
     }
     print OUT $s;
     close OUT;
     
     open OUT, '> '.$output_dir.'/'.$k."/info_seq.txt";
     $s=$chr.":".$sta."..".$end."\n".($#{$draw_pos{$k}}+1)."\n";
     for my $i(@{$draw_pos{$k}}){
          my @al=@{$drawl[$i]};my @ar=@{$drawr[$i]};
          my $qid=$qseqid[$al[0]];
          my $std="+";
          my $sseq=$h_seqs{$qid};
          if ($strd[$al[0]]==0) {
               $std="-";
               $sseq=~ tr/ATCG/TAGC/; $sseq=reverse $sseq;
          }      
          $s.=$qid."\n".$std."\n".$sseq."\n".($#al+1)."\t".($#ar+1)."\n";
          my $s1;
          my @l=@{$drawl_ad[$i]};my @r=@{$drawr_ad[$i]};
          for my $j(0..$#l){
               my($jl,$jr,$tem1,$tem2)=split(/_/,$l[$j]);
               $s1.=$jl."\t".($jr-$jl)."\n".get_misgap_p($al[$j])."\n";
               
          }
          for my $j(0..$#r){
               my($jl,$jr,$tem1,$tem2)=split(/_/,$r[$j]);
               $s1.=$jl."\t".($jr-$jl)."\n".get_misgap_p($ar[$j])."\n";
          }
          $s.=$s1;
     }
     print OUT $s;
     close OUT;
     
     draw_R($output_dir."/".$k.'/');
     for my $i(@{$draw_pos{$k}}){
          my @al=@{$drawl[$i]};my @ar=@{$drawr[$i]};
          my $qid=$qseqid[$al[0]];
          my $bits=0;
          for my $i(@al){$bits+=$bitscore[$i];}
           for my $i(@ar){$bits+=$bitscore[$i];}
          #$bits/=($#al+$#ar+1);
          push(@{$query_type{$qid}},"circle_pos,".$k.",".$bits);
      }  
}

sub print_gene_struct
{
     
   my ($s)= @_;
   open OUT_INFO, '> '.$output_dir.'/'.$s.'/info.txt';
   my ($chr,$sta,$end)=split(/_/,$s);
   my @choose_gene;my @choose_trans;
   for my $i(@{$chrom_gene{$chr}}){
       if($#{$gene_mRNA{$i}}==-1){next;}
       my $gsta=$gene_up{$i};
       my $gend=$gene_down{$i};
       #if ($i eq "LOC_Os05g37170") {my @aaa=@{$m_trans{"LOC_Os05g37170.5"}};$gsta=$gene_up{$i};}
       if (($sta>=$gsta && $sta<$gend) || ($end>$gsta && $end<=$gend)) {
             push(@choose_gene,$i);
             for my $j(@{$gene_mRNA{$i}}){
                  print_gene_struct_mRNA($s,$j,$sta,$end);
                  push(@choose_trans,$j);
               }
          }
     }
   if ($#choose_gene>=0) {
     print OUT_INFO ($#choose_trans+1)."\n".join("\t",@choose_trans)."\n".$s."\n",$gene_strd{$choose_gene[0]}."\n";
    # print OUT_INFO "1";
   }
   else{
     print OUT_INFO "1\nno_gene\n".$s."\n";
     open OUT2, '> '.$output_dir.'/'.$s."/no_gene";
     print OUT2 "intergenic\t".($end-$sta)."\n";
     close OUT2;
   }
   close OUT_INFO;
}

sub print_gene_struct_mRNA
{
     my ($kdir,$m,$sta,$end)=@_;
     open OUT, '>'.$output_dir.'/'.$kdir.'/'.$m;
     my %sort_trans;
     for my $i(@{$m_trans{$m}}){
         my($type,$msta,$mend)=split(/-/,$i);
         $sort_trans{$msta}=$i;
     }
     my $s="";my @a_type;my @a_sta;my @a_end;
    
     my $pre_end=0;
     for my $k(sort{$a<=>$b} keys %sort_trans){
          my($type,$msta,$mend)=split(/-/,$sort_trans{$k});
          
          if ($pre_end!=0) {
              if ($msta-$pre_end>0) {
                 push(@a_type,"intron");
                 push(@a_sta,$pre_end);
                 push(@a_end,$msta);
              }
          }
          push(@a_type,$type);
          push(@a_sta,$msta);
          push(@a_end,$mend);
          $pre_end=$mend;
     }
     if ($sta<$a_sta[0]) {$s.="intergenic\t".($a_sta[0]-$sta)."\n";}
     for my $i(0..$#a_type){
          if ($sta<=$a_sta[$i] && $end>=$a_end[$i]) {$s.=$a_type[$i]."\t".($a_end[$i]-$a_sta[$i])."\n";next;}
          if ($sta>=$a_sta[$i] && $end<=$a_end[$i]) {$s.=$a_type[$i]."\t".($end-$sta)."\n";next;}
          if ($sta<=$a_sta[$i] && $end<=$a_end[$i] && $end>$a_sta[$i]) {$s.=$a_type[$i]."\t".($end-$a_sta[$i])."\n";next;}
          if ($sta>=$a_sta[$i] && $sta<$a_end[$i] && $end>=$a_end[$i]) {$s.=$a_type[$i]."\t".($a_end[$i]-$sta)."\n";next;}
     }
     if ($end>$a_end[-1]) {$s.="intergenic\t".($end-$a_end[-1])."\n";}
     print OUT $s;
     close OUT;
     
}

sub print_gene_beside
{
     
     my ($s)=@_;
     open OUT_INFO, '>'.$output_dir.'/'.$s.'/infog.txt';
     my ($chr,$sta,$end)=split(/_/,$s);
     my @mid_gene;
     my $up_gene="";my $down_gene="";
     my $up_dis=9999999999;my $down_dis=9999999999;
     for my $i(@{$chrom_gene{$chr}}){
          my $gsta=$gene_up{$i};
          my $gend=$gene_down{$i};
          #if ($i eq "LOC_Os05g37170") {my @aaa=@{$m_trans{"LOC_Os05g37170.5"}};$gsta=$gene_up{$i};}
          if(($sta>=$gsta && $sta<$gend) || ($end>$gsta && $end<=$gend)) {
              push(@mid_gene,$i);
          }
          if ($gend<=$sta && $sta-$gend<$up_dis) {
              $up_dis=$sta-$gend;$up_gene=$i;
          }
          if ($gsta>=$end && $gsta-$end<$down_dis) {
              $down_dis=$gsta-$end;$down_gene=$i;
          }
     }
     my $l=9999999999;my $r=0;my $count=0;
     my $ml=9999999999;my $mr=0;
     for my $i(@mid_gene){
          $count++;
          if ($ml>$gene_up{$i}) {$ml=$gene_up{$i};}
          if ($mr<$gene_down{$i}) {$mr=$gene_down{$i};}  
     }
     if ($up_gene eq "") {
         if ($#mid_gene==-1) {$l=$gene_up{$down_gene};}
         else{$l=$ml;}
     }
     else{
         $l=$gene_up{$up_gene};$count++;
     }
      if ($down_gene eq "") {
         if ($#mid_gene==-1) {$r=$gene_down{$up_gene};}
         else{$r=$mr;}
     }
     else{
         $r=$gene_down{$down_gene};$count++;
     }
     
     print OUT_INFO $chr."\t".$sta."\t".$end."\n";
     my $out_info="";my $max_num=0;
     if (not $up_gene eq "") {
        my ($up_info,$num)=print_gene_beside_mRNA($up_gene);
        if ($num>$max_num) {$max_num=$num;} 
        $out_info.=$up_info;
     }
     for my $i(@mid_gene){
         my ($mid_info,$num)=print_gene_beside_mRNA($i);
         if ($num>$max_num) {$max_num=$num;} 
         $out_info.=$mid_info;
     }
     if (not $down_gene eq "") {
        my ($down_info,$num)=print_gene_beside_mRNA($down_gene);
        if ($num>$max_num) {$max_num=$num;} 
        $out_info.=$down_info;
     }
     print OUT_INFO $l."\t".$r."\n".$count."\t".$max_num."\n";
     print OUT_INFO $out_info;
     close OUT_INFO;
}

sub print_gene_beside_mRNA
{
     my ($sg)=@_;
     my $count=($#{$gene_mRNA{$sg}}+1);
     my $s_out=$sg."\t".$count."\t".$gene_strd{$sg}."\n";
     for my $i(@{$gene_mRNA{$sg}}){
          my %sort_trans;
          for my $j(@{$m_trans{$i}}){
              my($type,$msta,$mend)=split(/-/,$j);
              $sort_trans{$msta}=$j;
          }
          my @a_type;my $pre_end=0;
          for my $k(sort{$a<=>$b} keys %sort_trans){
              my($type,$msta,$mend)=split(/-/,$sort_trans{$k});  
              if ($pre_end!=0) {
                  if ($msta-$pre_end>0) {push(@a_type,"intron-".$pre_end."-".$msta);}
              }
              push(@a_type,$type."-".$msta."-".$mend);
              $pre_end=$mend;
          }
          my $st_out=$i."\t".($#a_type+1)."\n".join("\t",@a_type)."\n";
          $s_out.=$st_out;
     }
     return ($s_out,$count);
}

open OUT2, '> '.$output_dir.'/result.txt';
for my $k(keys %h_seq){
     print OUT2 $k."\n";
     if (exists($query_type{$k})) {
          my %ks;
          for my $i(@{$query_type{$k}}){
               my($type,$pos,$bits)=split(/,/,$i);
               $ks{$bits}=$type."\t".$pos."\t".$bits;
               #print OUT2 $k."\n";
          }
          for my $k(sort{$b<=>$a} keys %ks){
               print OUT2 $ks{$k}."\n";
          }
     }
     else{
          print OUT2 "no_hit\n";
     }
}

close OUT2;
sub get_misgap_p
{
     my($i)=@_;
     my $std=$strd[$i];
     my $qsq=uc($qseq[$i]);my $ssq=uc($sseq[$i]);
     my $mis=$mismatch[$i];my $gap=$gapopen[$i];
     my $qsta=$qstart[$i];
     my @cmis;my @cgap;
     if ($std==0) {
         $qsq=~ tr/ATCG/TAGC/;$qsq=reverse $qsq;
         $ssq=~ tr/ATCG/TAGC/;$ssq=reverse $ssq;
     }
     if ($mis>0 || $gap>0) {
         my $count=0;my $tagq=0;my $tags=0;
         for my $i(0..(length($ssq)-1)){
           if ((not substr($ssq,$i,1) eq substr($qsq,$i,1)) &&(not substr($qsq,$i,1) eq "-") && (not substr($ssq,$i,1) eq "-")){
               push(@cmis,($qsta+$i)."\t1");
           }
           if (substr($qsq,$i,1) eq "-" && $tagq==0){push(@cgap,($qsta+$i)."\t1");$tagq=1;}
           else{$tagq=0;}
           if (substr($ssq,$i,1) eq "-"&& $tags==0){push(@cgap,($qsta+$i)."\t1");$tags=1;}
           else{$tags=0;}
         } 
     }
     my $smis=$#cmis+1;
     for my $i(@cmis){$smis.="\n".$i;}
     my $sgap=$#cgap+1;
     for my $i(@cgap){$sgap.="\n".$i;}
     return $smis."\n".$sgap;
}
sub get_misgap_s
{
     my($i,$sta,$end)=@_;
     my $len=$end-$sta;
     my $std=$strd[$i];
     my $qsq=uc($qseq[$i]);my $ssq=uc($sseq[$i]);
     my $mis=$mismatch[$i];my $gap=$gapopen[$i];
     my $ssta=$sstart[$i];
     my @cmis;my @cgap;
     if ($std==0) {
         $qsq=~ tr/ATCG/TAGC/;$qsq=reverse $qsq;
         $ssq=~ tr/ATCG/TAGC/;$ssq=reverse $ssq;
     }
     if ($mis>0 || $gap>0) {
          my $count=0;
         for my $i(0..(length($ssq)-1)){
           if ((not substr($ssq,$i,1) eq substr($qsq,$i,1)) &&(not substr($qsq,$i,1) eq "-") && (not substr($ssq,$i,1) eq "-")){
               push(@cmis,($ssta+$i-$sta)."\t1\t".($len-($ssta+$i-$sta+1)));
           }
           if (substr($qsq,$i,1) eq "-" || substr($ssq,$i,1) eq "-"){
               push(@cgap,($ssta+$i-$sta)."\t1\t".($len-($ssta+$i-$sta+1)));
           }
         } 
     }
     my $smis=$#cmis+1;
     for my $i(@cmis){$smis.="\n".$i;}
     my $sgap=$#cgap+1;
     for my $i(@cgap){$sgap.="\n".$i;}
     return $smis."\n".$sgap;
}

sub getsum{
     my @a=@_;
     my $seql=$h_seq{$qseqid[$a[0]]};
     my %ht;
     for my $i(@a){
          for my $j($qstart[$i]..$qend[$i]){
               $ht{$j}=1;
          }
     }
     my $sum=0;
     for my $i(1..$seql){
          if (exists($ht{$i})) {
               $sum++;
          }    
     }
     return $sum/$seql;
}
sub get_bound_ad_l
{
     my ($tempa,$minu)=@_;
     my @a=@$tempa;
     my @qs;my @qe;my @ss;my @se;
     my $chr=$sseqid[$a[0]];
     for my $i(@a){
          push(@qs,$qstart[$i]);
          push(@qe,$qend[$i]);
          push(@ss,$sstart[$i]);
          push(@se,$send[$i]);
     }
     $qe[-1]-=$minu;
     $se[-1]-=$minu;
     for my $i(1..$#a){
          my $mi=$qe[$i-1]-$qs[$i];
          if ($mi<=0) {next;}
          for my $j(1..$mi){
               my $b=0;
               my $n1=$se[$i-1]-$j;
               my $n2=$ss[$i]+($mi-$j);
               if (exists($h_anno_up{$chr}{$n1}) || exists($h_anno_down{$chr}{$n1}) ) {$b++;}
               if (exists($h_anno_up{$chr}{$n2}) || exists($h_anno_down{$chr}{$n2}) ) {$b++;}
               if ($b==2) {
                    $se[$i-1]-=$j;$qe[$i-1]-=$j;
                    $ss[$i]+=($mi-$j);$qs[$i]+=($mi-$j);
               }  
          }
     }
     my @aa;
     for my $i(0..$#a){
          push(@aa,$qs[$i]."_".$qe[$i]."_".$ss[$i]."_".$se[$i]);
     }
     return @aa;
}

sub get_bound_ad_r
{
     my @a=@_;
     my @qs;my @qe;my @ss;my @se;
     my $chr=$sseqid[$a[0]];
     for my $i(@a){
          push(@qs,$qstart[$i]);
          push(@qe,$qend[$i]);
          push(@ss,$sstart[$i]);
          push(@se,$send[$i]);
     }
     for my $i(1..$#a){
          my $mi=$qe[$i]-$qs[$i-1];
          if ($mi<=0) {next;}
          for my $j(1..$mi){
               my $b=0;
               my $n1=$se[$i]-$j;
               my $n2=$ss[$i-1]+($mi-$j);
               if (exists($h_anno_up{$chr}{$n1}) || exists($h_anno_down{$chr}{$n1}) ) {$b++;}
               if (exists($h_anno_up{$chr}{$n2}) || exists($h_anno_down{$chr}{$n2}) ) {$b++;}
               if ($b==2) {
                    $se[$i]-=$j;$qe[$i]-=$j;
                    $ss[$i-1]+=($mi-$j);$qs[$i-1]+=($mi-$j);
               }  
          }
     }
      my @aa;
     for my $i(0..$#a){
          push(@aa,$qs[$i]."_".$qe[$i]."_".$ss[$i]."_".$se[$i]);
     }
     return @aa;
}

sub check_splice_site_two
{
     my($sta,$end,$chr)=@_;
     my $b=0;
     if (exists($h_anno_up{$chr}{$sta}) || exists($h_anno_down{$chr}{$sta}) ){$b++};
     if (exists($h_anno_up{$chr}{$end}) || exists($h_anno_down{$chr}{$end}) ){$b++};
     if ($b==2) {return 0;}
     $h_refer{$chr}=uc($h_refer{$chr});
     my $sd=substr($h_refer{$chr},$sta-2,2);
     my $ed=substr($h_refer{$chr},$end,2);
     if (($sd eq "AG" && $ed eq "GT") || ($sd eq "AC" && $ed eq "CT")) {
         return 0;
     }
     
     my $look_for=20;
     for my $i(1..20){
          my $sseq=substr($h_refer{$chr},$sta,$i);
          my $eseq=substr($h_refer{$chr},$end,$i);
          if (not $sseq eq $eseq ) {last;} 
          my $sd=substr($h_refer{$chr},$sta+$i-2,2);
          my $ed=substr($h_refer{$chr},$end+$i,2);
           if (($sd eq "AG" && $ed eq "GT") || ($sd eq "AC" && $ed eq "CT")) { return $i;}
     }
     for my $j(1..20){
          my $sseq=substr($h_refer{$chr},$sta-$j,$j);
          my $eseq=substr($h_refer{$chr},$end-$j,$j);
          if (not $sseq eq $eseq ) {last;} 
          my $sd=substr($h_refer{$chr},$sta-$j-2,2);
          my $ed=substr($h_refer{$chr},$end-$j,2);
           if (($sd eq "AG" && $ed eq "GT") || ($sd eq "AC" && $ed eq "CT")) { return -$j;}
     }
     return 0;
}
sub check_splice_site
{
     my ($aa1,$aa2,$chr)=@_;
   #my (@a)=@_;
   #my @a1=$a[0];
  # my @a2=$a[1];
   #my $chr=$a[2];
  # my $ss1=$a1[-1];
  # my $ss2=$a2[-1];
  my @a11=@$aa1;my @a22=@$aa2;
   my @c1=split(/_/,$a11[-1]);
   my @c2=split(/_/,$a22[-1]);
   my $move=check_splice_site_two($c2[2],$c1[3],$chr);
   if ($move!=0) {
     #$a11[-1]=$c2[0]."_".$c2[1]."_".($c2[2]+$move)."_".$c2[3];
     #$a22[-1]=$c2[0]."_".$c2[1]."_".$c2[2]."_".($c2[3]+$move);
     $a11[-1]=$c1[0]."_".($c1[1]+$move)."_".$c1[2]."_".($c1[3]+$move);
     $a22[-1]=($c2[0]+$move)."_".$c2[1]."_".($c2[2]+$move)."_".$c2[3];
   }
   my $k=$chr."_".($c2[2]+$move)."_".($c1[3]+$move);
   push(@{$draw_pos{$k}},$count_draw);
   $count_draw++;
   return (\@a11,\@a22);
}

sub draw_R
{
    my ($dir)=@_;
    my $dir_f=$my_dir."/format.txt";
    my $dir_r=$my_dir."/R_script";
    if($color_type eq "3D"){
         $command=qq(Rscript  $dir_r/draw_seq.R  $dir $dir_f  $dir/Seq.pdf);
         system($command);
         $command=qq(Rscript  $dir_r/3D_circ_seq.R  $dir $dir_f  $dir/circle_seq.pdf);
         system($command);
         $command=qq(Rscript  $dir_r/3D_linear_seq.R  $dir  $dir_f  $dir/linear_seq.pdf);
         system($command);
      }
    else{
        $command=qq(Rscript  $dir_r/draw_seq.R  $dir $dir_f  $dir/Seq.pdf);
         system($command);
         $command=qq(Rscript  $dir_r/2D_circ_seq.R  $dir $dir_f  $dir/circle_seq.pdf);
         system($command);
         $command=qq(Rscript  $dir_r/2D_linear_seq.R  $dir  $dir_f  $dir/linear_seq.pdf);
         system($command);
    }
}

sub draw_linear_R
{
    my ($dir)=@_;
    my $dir_f=$my_dir."/format.txt";
    my $dir_r=$my_dir."/R_script";
    if($color_type eq "3D"){
         $command=qq(Rscript  $dir_r/3D_linear_seq.R  $dir  $dir_f  $dir/linear_seq.pdf);
         system($command);
      }
    else{
         $command=qq(Rscript  $dir_r/2D_linear_seq.R  $dir  $dir_f  $dir/linear_seq.pdf);
         system($command);
    }
}

print "\n";