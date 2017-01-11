#!/usr/bin/perl -w
use strict;
use File::Basename;
my $my_dir=dirname(__FILE__);
my $input_config =$my_dir."/pos_config.txt";



my %h_anno_up;
my %h_anno_down;


my %chrom_gene;
my %gene_mRNA;
my %gene_up;
my %gene_down;
my %gene_strd;
my %m_trans;

my $input_anno;
my $anno_type;
my $input_chrom;
my $input_sta;
my $input_end;

my $color_type="3D";
my $output_dir="circRNA_pocket_pos_output";

my $b_test_config=test_config();
my $b_test_format=test_format();
if($b_test_config!=1){die "Exit.";}
if($b_test_format!=1){die "Exit.";}
sub test_config
{
    open IN_config, '<'.$input_config or die "Can't open $input_config. $!";
    my @c_name=("query_chrom","query_start_pos","query_end_pos","annotation_file","annotation_type","color_type","output_dir");
    my @b_name=(0,0,0,0,0,0,0);
    while (<IN_config>){
          my $s=$_;
          chomp($s);
          my @c=split(/=/,$s);
          if (substr($s,0,1) eq "#" || $s eq "") {next;}
     
          if($c[0] eq "query_chrom"){
             $b_name[0]=1;
             $input_chrom=$c[1];next;
          }
          if($c[0] eq "query_start_pos"){
             $b_name[1]=1;
             $input_sta=int($c[1]);next;
          }
          if($c[0] eq "query_end_pos"){
             $b_name[2]=1;
             $input_end=int($c[1]);next;
          }
          if($c[0] eq "annotation_file"){
             $b_name[3]=1;
             $input_anno = $c[1];next;
          }
          if($c[0] eq "annotation_type"){
             $b_name[4]=1;
             $anno_type = $c[1];next;
          }
          if($c[0] eq "color_type"){
             $b_name[5]=1;
             $color_type = $c[1];next;
          }
          if($c[0] eq "output_dir"){
             $b_name[6]=1;
             $output_dir = $c[1];next;
          }
    }
    close IN_config;
    for my $i(0..4){
         if($b_name[$i]==0){
            print "$input_config. Can not find the item $c_name[$i].\n";
            return 0;
         }
    }
    if($input_sta<0){
          print "$input_config. query_start_pos should be a positive integer.\n";
          return 0;
    }
    if($input_end<0){
          print "$input_config. query_end_pos should be a positive integer.\n";
          return 0;
    }
    if($input_sta>$input_end){
          print "$input_config. query_end_pos should be bigger than query_start_pos.\n";
          return 0;
    }
    if($input_end-$input_sta>2000000){
          print "$input_config. query_end_pos - query_start_pos should be no more than 2000000.\n";
          return 0;
    }

    return 1;
}
sub test_format
{

     my $input_format =$my_dir."/pos_format.txt";
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
              my @c=split(/\t/,$s);
              my @c1=split(/_/,$c[0]);my @c2=split(/_/,$c[1]);
              $x1=$c1[0];$y1=$c1[1];$x2=$c2[0];$y2=$c2[1];
              if($x1<2 || $x2<1 || $x1>10 || $x2>10 || $y2<=0 || $y2>10 || $y1<2 || $y1>6){
                  print "$input_format\tline 2.x1_y1 x2_y2 should be (2 <= x1 <= 10; 1 <= x2 <= 10;0 < y2 <= 10;2 <= y1 <=6)\n";
                  return 0;
              }
              next;
          }
          if($count==3){
              if((not $s eq "unclockwise") && (not $s eq "clockwise")){
                  print "$input_format\tline 3.It should be (unclockwise;clockwise)\n";return 0;
              }
              next;
          }
          if($count==4){
              if((not $s eq "show_UTR") && (not $s eq "no_UTR")){
                  print "$input_format\tline 4.It should be (show_UTR;no_UTR)\n";return 0;
              }
              next;
          }
           if($count==5){
              if((not $s eq "show_gene_direction") && (not $s eq "no_gene_direction")){
                  print "$input_format\tline 5.It should be (show_gene_direction;no_gene_direction)\n";return 0;
              }
              next;
          }

     }
     close IN_format;
     return 1;
}


mkdir($output_dir);
open IN_anno, '<'.$input_anno or die "Can't open $input_anno. $!";

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
    if ($type eq "exon") {
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
    if ($type eq "mRNA" || $type eq "transcript") {
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
          $r=~s/:/_/;
          if ($l eq $gt) {return $r; }       
     }
     return;
  }

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



make_info_file($input_chrom."_".$input_sta."_".$input_end);

sub make_info_file
{
     my ($k)=@_;
     mkdir($output_dir."/".$k);
    # print_gene_struct($k);
     print_gene_struct($k);
     print_gene_beside($k);
     draw_R($output_dir."/".$k.'/');
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
      # print "myprint0 i:".$i."\tgsta:".$gsta."\tgend:".$gend."\n";
       #if ($i eq "LOC_Os05g37170") {my @aaa=@{$m_trans{"LOC_Os05g37170.5"}};$gsta=$gene_up{$i};}
       if (($sta>=$gsta && $sta<$gend) || ($end>$gsta && $end<=$gend)) {
             push(@choose_gene,$i);
             for my $j(@{$gene_mRNA{$i}}){
                  #print "myprint1 i:".$i."\tj:".$j."\ts:".$s."\tsta:".$sta."\tend:".$end."\n";
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
     #print "myprint2 kdir:".$kdir."\tm:".$m."\tsta:".$sta."\tend:".$end."\n";
     open OUT, '>'.$output_dir.'/'.$kdir.'/'.$m;
     #print "myprint:".$output_dir.'/'.$kdir.'/'.$m."\n";
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
          if((not exists($gene_mRNA{$i})) || $#{$gene_mRNA{$i}}==-1){next;}
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
         if($ml<$l){$l=$ml-5;}
     }
      if ($down_gene eq "") {
         if ($#mid_gene==-1) {$r=$gene_down{$up_gene};}
         else{$r=$mr;}
     }
     else{
         $r=$gene_down{$down_gene};$count++;
         if($mr>$r){$r=$mr+5;}
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



sub draw_R
{
    my ($dir)=@_;
    my $dir_f=$my_dir."/pos_format.txt";
    my $dir_r=$my_dir."/R_script";
    my $command;
    if($color_type eq "3D"){
         $command=qq(Rscript  $dir_r/3D_circ.R  $dir $dir_f  $dir/circle.pdf);
         system($command);
         $command=qq(Rscript  $dir_r/3D_linear.R  $dir  $dir_f  $dir/linear.pdf);
         system($command);
      }
    else{
  
         $command=qq(Rscript  $dir_r/2D_circ.R  $dir $dir_f  $dir/circle.pdf);
         system($command);
         $command=qq(Rscript  $dir_r/2D_linear.R  $dir  $dir_f  $dir/linear.pdf);
         system($command);
    }
}


print "\n";