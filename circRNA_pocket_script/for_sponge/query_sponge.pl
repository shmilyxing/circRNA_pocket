#!/usr/bin/perl -w
use strict;
use File::Basename;
my $my_dir=dirname(__FILE__);
my $input_config =$my_dir."/sponge_config.txt";



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
my $input_genome;
my $input_microRNA;

my $software_type;
my $RNAhybrid_s="";
my $RNAhybrid_p="";
my $RNAhybrid_e="";

my $miranda_en="";
my $miranda_sc="";
my $color_type="3D";
my $output_dir="circRNA_pocket_sponge_output";

my $b_test_config=test_config();
my $b_test_format=test_format();
if($b_test_config!=1){die "Exit.";}
if($b_test_format!=1){die "Exit.";}
sub test_config
{
    open IN_config, '<'.$input_config or die "Can't open $input_config. $!";
    my @c_name=("query_chrom","query_start_pos","query_end_pos","annotation_file","annotation_type","gemone_reference","microRNA_sponge_detection_software","mature_microRNA_file");
    my @b_name=(0,0,0,0,0,0,0,0);
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
          if($c[0] eq "gemone_reference"){
             $b_name[5]=1;
             $input_genome = $c[1];next;
          }
          if($c[0] eq "microRNA_sponge_detection_software"){
             $b_name[6]=1;
             $software_type = $c[1];next;
          }
          if($c[0] eq "mature_microRNA_file"){
             $b_name[7]=1;
             $input_microRNA = $c[1];next;
          }
          if($c[0] eq "color_type"){$color_type = $c[1];next;}
          if($c[0] eq "output_dir"){$output_dir = $c[1];next;}

          if($c[0] eq "RNAhybrid-s"){$RNAhybrid_s = $c[1];next;}
          if($c[0] eq "RNAhybrid-p"){$RNAhybrid_p = $c[1];next;}
          if($c[0] eq "RNAhybrid-e"){$RNAhybrid_e = $c[1];next;}

          if($c[0] eq "miranda-en"){$miranda_en = $c[1];next;}
          if($c[0] eq "miranda-sc"){$miranda_sc = $c[1];next;}
    }
    close IN_config;
    for my $i(0..6){
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
    if((not $software_type eq "RNAhybrid") && (not $software_type eq "miranda")){
         print "$input_config. microRNA_sponge_detection_software - the software used to detect microRNA sponge should be (RNAhybrid;miranda).\n";
          return 0;
    }
    else{
           my $reg1 = qr/^-?\d+(\.\d+)?$/;
           my $reg2 = qr/^-?0(\d+)?$/;
          if($software_type eq "RNAhybrid"){
               if((not $RNAhybrid_s eq "3utr_fly") && (not $RNAhybrid_s eq "3utr_worm") && (not $RNAhybrid_s eq "3utr_human")){
                     print "$input_config. RNAhybrid_s should be (3utr_fly;3utr_worm;3utr_human)\n";return 0;
               }
               
               if($RNAhybrid_p =~ $reg1 && $RNAhybrid_p !~ $reg2 && $RNAhybrid_p>0 && $RNAhybrid_p<1){;}
               else{print "$input_config. RNAhybrid-p should be (0<p<1)\n";return 0; }
               if($RNAhybrid_e =~ $reg1 && $RNAhybrid_e !~ $reg2 && $RNAhybrid_e>0){;}
               else{print "$input_config. RNAhybrid-e should be (e>0)\n";return 0; }
           }
       
         else{
              if($software_type eq "miranda"){
                  if($miranda_en =~ $reg1 && $miranda_en !~ $reg2 && $miranda_en < 0){;}
                  else{print "$input_config. miranda-en should be(en<0)\n";return 0; }
                  if($miranda_sc =~ $reg1 && $miranda_sc !~ $reg2 && $miranda_sc > 0){;}
                  else{print "$input_config. miranda-sc should be (sc>0)\n";return 0; }
              }
          }
    }
    return 1;
}
sub test_format
{

     my $input_format =$my_dir."/sponge_format.txt";
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
               if($s eq "Blue" || $s eq "Green" || $s eq "Yellow" || $s eq "Red" || $s eq "Purple" || $s eq "Brown" || $s eq "Grey" || $s eq "Pink" || $s eq "Gold")
               {;}
               else{
                  print "$input_format\tline 1.The color of microRNA sponge should be (Blue;Green;Yellow;Red;Purple;Brown;Grey;Pink;Gold)\n";
                  return 0;
               }
               next;
           }
          if($count==3){
              my @c=split(/\t/,$s);
              my @c1=split(/_/,$c[0]);my @c2=split(/_/,$c[1]);
              $x1=$c1[0];$y1=$c1[1];$x2=$c2[0];$y2=$c2[1];
              if($x1<2 || $x2<1 || $x1>10 || $x2>10 || $y2<=0 || $y2>10 || $y1<2 || $y1>6){
                  print "$input_format\tline 2.x1_y1 x2_y2 should be (2 <= x1 <= 10; 1 <= x2 <= 10;0 < y2 <= 10;2 <= y1 <=6)\n";
                  return 0;
              }
              next;
          }
          if($count==4){
              if((not $s eq "unclockwise") && (not $s eq "clockwise")){
                  print "$input_format\tline 3.It should be (unclockwise;clockwise)\n";return 0;
              }
              next;
          }
          if($count==5){
              if((not $s eq "show_UTR") && (not $s eq "no_UTR")){
                  print "$input_format\tline 4.It should be (show_UTR;no_UTR)\n";return 0;
              }
              next;
          }
           if($count==6){
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
open IN_refer, '<'.$input_genome or die "Can't open $input_genome. $!";
print "Reading genome annotation file..\n";
if ($anno_type eq "gff3") {get_anno_gff3();}
if ($anno_type eq "gtf") {get_anno_gtf();}
if ($anno_type eq "refFlat") {get_anno_refFlat();}


print "Reading genome file:\n";
my $s4;
my %h_refer;
my $b_get_chrom=0;
while (<IN_refer>){
      my $s=$_;
     chomp($s);
     my $s3=substr($s,0,1);
    if ($s3 eq ">") {
       my @c=split(/ /,$s);
       $s4=substr($c[0],1);
       if($b_get_chrom==1){last;}
       if($input_chrom eq $s4){$b_get_chrom=1;}
       print $s4."\n";
    }
    else
    {    
       $h_refer{$s4}.=$s;      
    }   
}
close IN_refer;

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
     my $dir=$output_dir."/".$k;
     my $command;my $b_command_res;
     if($software_type eq "RNAhybrid"){
          $command=qq(RNAhybrid -c -e -$RNAhybrid_e -p $RNAhybrid_p -s $RNAhybrid_s -t $dir/for_micro.fa -q $input_microRNA > $dir/RNAhybrid_res.out);
          $b_command_res=system($command);
          if($b_command_res!=0){die "RNAhybrid  failed!";}
          print_micro_struct_RNAhybrid($k);
     }
     if($software_type eq "miranda"){
          $command=qq(miranda  $input_microRNA $dir/for_micro.fa  -en $miranda_en -sc $miranda_sc -out $dir/miranda.out -keyval -quiet);
          $b_command_res=system($command);
          if($b_command_res!=0){die "miranda  failed!";}
          $command=qq(grep '>' miranda.out > miranda_res.out);
          $b_command_res=system($command);
          if($b_command_res!=0){die "miranda grep failed!";}
          print_micro_struct_miranda($k);
     }
     
     draw_R($output_dir."/".$k.'/');
}


my %h_ei_trans;
sub print_micro_struct_RNAhybrid
{
     
      my ($k)= @_;
      open OUT_MICRO, '> '.$output_dir.'/'.$k.'/info_micro.txt';
      my ($chr,$sta,$end)=split(/_/,$k);
      my %h_micro_struct;
      my $input_micro=$output_dir.'/'.$k.'/RNAhybrid_res.out';
      open IN_micro, '<'.$input_micro or die "Can't open $input_micro. $!";
      my $out_s="";
      while (<IN_micro>){
           my $s=$_;
           chomp($s);
           my @c = split(/:/,$s);
           my @c2 = split(/;/,$c[0]);
           my $len = $c[3];my $pos= $c[6];my $all_len=0;
           my $trans_id=$c2[1];my $micro_id=$c[2];
           if($trans_id eq "all"){
                #$out_s.=$trans_id.";".$micro_id.";".$pos."\t".$len."\t".($end-$sta-$len-$pos)."\n" ;
                 if(exists($h_micro_struct{$trans_id})){
                        $h_micro_struct{$trans_id}.=";".$pos."\t".$len."\t".($end-$sta-$len-$pos)."\t".$micro_id;
                  }else{
                        $h_micro_struct{$trans_id}.=$pos."\t".$len."\t".($end-$sta-$len-$pos)."\t".$micro_id;
                  }
                next;
           }
           my ($s_intron,$s_exon)=split(/\t/,$h_ei_trans{$trans_id});
           my @c_intron=split(/_/,$s_intron);my @c_exon=split(/_/,$s_exon);
           for my $i(0..$#c_exon){
              if($pos<$c_exon[$i] && $i==0){last;}
              if($pos<$c_exon[$i] && $i>0 && $pos>=$c_exon[$i-1]){
                $pos+=$c_intron[$i-1];
                last;
              }
           }
            if(exists($h_micro_struct{$trans_id})){
                       $h_micro_struct{$trans_id}.=";".$pos."\t".$len."\t".($end-$sta-$len-$pos)."\t".$micro_id;
            }else{
                       $h_micro_struct{$trans_id}.=$pos."\t".$len."\t".($end-$sta-$len-$pos)."\t".$micro_id;
            }
           #$out_s.=$trans_id.";".$micro_id.";".$pos."\t".$len."\t".($end-$sta-$len-$pos)."\n" ;
      }
      close IN_micro;
      for my $k(keys %h_micro_struct){
            my @c=split(/;/,$h_micro_struct{$k});
            print OUT_MICRO $k.";".($#c+1).";".$h_micro_struct{$k}."\n";
      }
      #print OUT_MICRO $out_s;
      print OUT_MICRO "---\n";
      for my $k(keys %h_micro_struct){
              my @c=split(/;/,$h_micro_struct{$k});
              #@c=splice(@c,2);
              my $trans_s=print_micro_struct_get_info(join(";",@c));   
              print OUT_MICRO $k.":\n".$trans_s."\n";
      }
      close OUT_MICRO;
}

sub print_micro_struct_get_info
{
     
   my ($s)= @_;
   my @c=split(/;/,$s);
   my %h1;my %h2;

   for my $i(@c){
     my @c2=split(/\t/,$i);
     $h1{$c2[3]}=$c2[0];
     $h2{$c2[3]}=$c2[1]+$c2[0];
   }
   my $pre=-1;my $ss="";
   for my $k(sort{$h1{$a}<=>$h1{$b}} keys %h1){

        if($pre==-1){$ss=$k;$pre=$h2{$k};next;}
        if($h1{$k}<=$pre){$ss.=" ".$k;$pre=$h2{$k};}
        else{
            $ss.="\n".$k;$pre=$h2{$k};
        }
   }
   return $ss."\n";
}


sub print_micro_struct_miranda
{
     
      my ($k)= @_;
      open OUT_MICRO, '> '.$output_dir.'/'.$k.'/info_micro.txt';
      my ($chr,$sta,$end)=split(/_/,$k);
      my %h_micro_struct;
      my $input_micro=$output_dir.'/'.$k.'/miranda_res.out';
      open IN_micro, '<'.$input_micro or die "Can't open $input_micro. $!";
      my $out_s="";
      while (<IN_micro>){
           my $s=$_;
           chomp($s);
           my @c = split(/\t/,$s);
           my $len = $c[7];my $pos= $c[9];my $all_len=0;
           my $trans_id=$c[1];my $micro_id=substr($c[0],2);
           my @c_pos= split(/ /,$pos);
           if($trans_id eq "all"){
                #$out_s.=$trans_id.";".$micro_id.";".$pos."\t".$len."\t".($end-$sta-$len-$pos)."\n" ;
                for my $i(@c_pos){
                    if(exists($h_micro_struct{$trans_id})){
                        $h_micro_struct{$trans_id}.=";".$i."\t".$len."\t".($end-$sta-$len-$i)."\t".$micro_id;
                    }else{
                        $h_micro_struct{$trans_id}.=$i."\t".$len."\t".($end-$sta-$len-$i)."\t".$micro_id;
                    }
                }
                next;
           }
           my ($s_intron,$s_exon)=split(/\t/,$h_ei_trans{$trans_id});
           my @c_intron=split(/_/,$s_intron);my @c_exon=split(/_/,$s_exon);
           for my $i(0..$#c_exon){
                for my $i(0..$#c_pos){
                    if($c_pos[$i] < $c_exon[$i] && $i==0){last;}
                    if($c_pos[$i] < $c_exon[$i] && $i>0 && $c_pos[$i] >= $c_exon[$i-1]){
                         $c_pos[$i]+=$c_intron[$i-1];last;
                    }
                }
           }
           for my $i(@c_pos){
                if(exists($h_micro_struct{$trans_id})){
                       $h_micro_struct{$trans_id}.=";".$i."\t".$len."\t".($end-$sta-$len-$i)."\t".$micro_id;
                 }else{
                       $h_micro_struct{$trans_id}.=$i."\t".$len."\t".($end-$sta-$len-$i)."\t".$micro_id;
                 }
           }
           
           #$out_s.=$trans_id.";".$micro_id.";".$pos."\t".$len."\t".($end-$sta-$len-$pos)."\n" ;
      }
      close IN_micro;
      for my $k(keys %h_micro_struct){
            my @c=split(/;/,$h_micro_struct{$k});
            print OUT_MICRO $k.";".($#c+1).";".$h_micro_struct{$k}."\n";
      }
      #print OUT_MICRO $out_s;
      print OUT_MICRO "---\n";
      for my $k(keys %h_micro_struct){
              my @c=split(/;/,$h_micro_struct{$k});
              #@c=splice(@c,2);
              my $trans_s=print_micro_struct_get_info(join(";",@c));   
              print OUT_MICRO $k.":\n".$trans_s."\n";
      }
      close OUT_MICRO;
}


sub print_gene_struct
{
     
   my ($s)= @_;
   open OUT_INFO, '> '.$output_dir.'/'.$s.'/info.txt';
   open OUT_SEQ, '> '.$output_dir.'/'.$s.'/for_micro.fa';
   my ($chr,$sta,$end)=split(/_/,$s);
   my @choose_gene;my @choose_trans;
   my $circ_seq="";
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
                  $circ_seq.=print_gene_struct_seq($s,$j,$sta,$end);
                  $h_ei_trans{$j}=print_gene_struct_ei($s,$j,$sta,$end);
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
   if($circ_seq eq ""){
      $circ_seq=">".$s.";all\n".substr($h_refer{$chr},$sta,$end-$sta)."\n";
   }
   print OUT_SEQ $circ_seq;
   close OUT_SEQ;
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

sub print_gene_struct_seq
{
     my ($kdir,$m,$sta,$end)=@_;
     my ($chr,$ksta,$kend)=split(/_/,$kdir);
     my %sort_trans;
     for my $i(@{$m_trans{$m}}){
         my($type,$msta,$mend)=split(/-/,$i);
         $sort_trans{$msta}=$i;
     }
     my $s="";my @a_type;my @a_sta;my @a_end;
    
     my $pre_end=0;
     for my $k(sort{$a<=>$b} keys %sort_trans){
          my($type,$msta,$mend)=split(/-/,$sort_trans{$k});
          push(@a_type,$type);
          push(@a_sta,$msta);
          push(@a_end,$mend);
     }
     my $b_sta=0;my $b_end=0;my $s_seq=">".$kdir.";".$m."\n";
     for my $i(0..$#a_type){
          if($a_sta[$i] == $sta){$b_sta=1;}
          if($a_end[$i] == $end){$b_end=1;}
          if ($sta<=$a_sta[$i] && $end>=$a_end[$i]) {$s_seq.=substr($h_refer{$chr},$a_sta[$i],$a_end[$i]-$a_sta[$i]);next;}
          if ($sta>=$a_sta[$i] && $end<=$a_end[$i]) {$s_seq.=substr($h_refer{$chr},$sta,$end-$sta);next;}
          if ($sta<=$a_sta[$i] && $end<=$a_end[$i] && $end>$a_sta[$i]) {$s_seq.=substr($h_refer{$chr},$a_sta[$i],$end-$a_sta[$i]);next;}
          if ($sta>=$a_sta[$i] && $sta<$a_end[$i] && $end>=$a_end[$i]) {$s_seq.=substr($h_refer{$chr},$sta,$a_end[$i]-$sta);next;}
     }

    if($b_sta==0 || $b_end==0){
         return "";
    }
     else{
        return $s_seq."\n";
     }
}

sub print_gene_struct_ei
{
     my ($kdir,$m,$sta,$end)=@_;
     #print "myprint2 kdir:".$kdir."\tm:".$m."\tsta:".$sta."\tend:".$end."\n";
    
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
              if ($msta-$pre_end>=0) {
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
     my $b_sta=0;my $b_end=0;
     my @intron;my @exon;my $intron_all=0;my $exon_all=0;
     for my $i(0..$#a_type){
          if($a_sta[$i] == $sta){$b_sta=1;}
          if($a_end[$i] == $end){$b_end=1;}
          if ($sta<=$a_sta[$i] && $end>=$a_end[$i]) {
              if($a_type[$i] eq "intron"){$intron_all+=$a_end[$i]-$a_sta[$i];push(@intron,$intron_all);}
              else{$exon_all+=$a_end[$i]-$a_sta[$i];push(@exon,$exon_all);}
          }
          
     }
   
    if($b_sta==0 || $b_end==0){return "";}
     else{
        return join("_",@intron)."\t".join("_",@exon);
     }
     
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
    my $dir_f=$my_dir."/sponge_format.txt";
    my $dir_r=$my_dir."/R_script";
    my $command;
    if($color_type eq "3D"){
         $command=qq(Rscript  $dir_r/3D_circ_micro.R  $dir $dir_f  $dir/circ_micro.pdf);
         system($command);
         $command=qq(Rscript  $dir_r/3D_linear_micro.R  $dir  $dir_f  $dir/linear.pdf);
         system($command);
      }
    else{
  
         $command=qq(Rscript  $dir_r/2D_circ_micro.R  $dir $dir_f  $dir/circ_micro.pdf);
         system($command);
         $command=qq(Rscript  $dir_r/2D_linear_micro.R  $dir  $dir_f  $dir/linear.pdf);
         system($command);
    }
}


print "\n";