seq_entraper = 'GGCCTTCGGGCCAAGAAAGCCGCGGCTGTAAATAAAGGCAAATCCGCCACTCTCTGAGGGGCCGTTTGTCCTCCCCTCCTTGGCCCCGTCCTTCCCGCCGCCCCCTCCCGGCTCCCGGGCCCCGCGGCGCCCCGGCCCGAGCTCCTCCATTTAATCGGATTTGGGAGAAGGGGAGGATAAATCACGGCAGCAGCTTTACGGTCCCGGAGGAGAGGCGAGCCGCAGACAGGCACACCCCGGCCGGCGATAAAAACCGCCGCTGAAAGCCCACGGAGCAATTTCCCGGGACCCCGAGCGACGCCATTACAGGAATGTAATTTTGCCCGGATGAGGCCCCGAGTTTAATTATCCTCGCGGAGGAATTTCAATGCG';
[requestID ,requestTime] = blastncbi(seq_entraper,'blastn');
blast_data = getblast(requestID,'WaitTime',requestTime);
[start,~] = blast_data.Hits(4).HSPs.SubjectIndices; 
first_match = strsplit(blast_data.Hits(4).Name,'|');
accession = char(first_match(4));
entraper_data = getgenbank(accession);
seq_entraper = entraper_data.Sequence(start:start+790);
mi339 = 'CCGGCTCCCGGGCCCCGCGGCGCC';
mi663b = 'CAGGCACACCCCGGCCGGC';
mi95 = 'TCAACGGGTATTTATTGAGCA';
[score_1,align_1,start_1] = swalign(seq_entraper,mi339,'Alphabet','nt','Showscore',true);
[score_2,align_2,start_2] = swalign(seq_entraper,mi663b,'Alphabet','nt','Showscore',true);
[score_3,align_3,start_3] = swalign(seq_entraper,mi95,'Alphabet','nt','Showscore',true);
[score_0,align_0,start_0] = swalign(seq_entraper,cat(2,mi663b,mi339,mi95),'Alphabet','nt','Showscore',true);

%%The code for normal txt file which cause the program to lose data
% _filename1 = '5bp_deletion_1.txt';
% fid = fopen(filename1,'w');
% fprintf(fid,cat(2,seq_entraper(1),seq_entraper(7:length(seq_entraper)),';'));
% for ii = 2:82
    % new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)),';');
    % fid1 = fopen(filename1,'a');
    % fprintf(fid1,new_seq);
% end
% fclose all;
% filename2 = '5bp_deletion_2.txt';
% fid = fopen(filename2,'w');
% fprintf(fid,cat(2,seq_entraper(1:107),seq_entraper(113:length(seq_entraper)),';'));
% for ii = 108:136
    % new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)),';');
    % fid2 = fopen(filename2,'a');
    % fprintf(fid2,new_seq);
% end
% fclose all;
% filename3 = '5bp_deletion_3.txt';
% fid = fopen(filename3,'w');
% fprintf(fid,cat(2,seq_entraper(1:160),seq_entraper(166:length(seq_entraper)),';'));
% for ii = 161:300
    % new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)),';');
    % fid3 = fopen(filename3,'a');
    % fprintf(fid3,new_seq);
% end
% fclose all;
% filename4 = '5bp_deletion_4.txt';
% fid = fopen(filename4,'w');
% fprintf(fid,cat(2,seq_entraper(1:301),seq_entraper(307:length(seq_entraper)),';'));
% for ii = 302:452
    % new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)),';');
    % fid4 = fopen(filename4,'a');
    % fprintf(fid4,new_seq);
% end
% fclose all;
% filename5 = '5bp_deletion_5.txt';
% fid = fopen(filename5,'w');
% fprintf(fid,cat(2,seq_entraper(1:453),seq_entraper(459:length(seq_entraper)),';'));
% for ii = 454:532
    % new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)),';');
    % fid5 = fopen(filename5,'a');
    % fprintf(fid5,new_seq);
% end
% fclose all;
% filename6 = '5bp_deletion_6.txt';
% fid = fopen(filename6,'w');
% fprintf(fid,cat(2,seq_entraper(1:554),seq_entraper(560:length(seq_entraper)),';'));
% for ii = 555:685
    % new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)),';');
    % fid6 = fopen(filename6,'a');
    % fprintf(fid6,new_seq);
% end
% fclose all;
% filename7 = '5bp_deletion_7.txt';
% fid = fopen(filename7,'w');
% fprintf(fid,cat(2,seq_entraper(1:686),seq_entraper(692:length(seq_entraper)),';'));
% for ii = 687:length(seq_entraper)
    % new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)),';');
    % fid7 = fopen(filename7,'a');
    % fprintf(fid7,new_seq);
% end
% fclose all;_



filename1 = '4bp_deletion_1_1.txt';
for ii = 1:140
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+5:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename1,seqname,new_seq);
end
filename2 = '4bp_deletion_1_2.txt';
for ii = 141:280
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+5:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename2,seqname,new_seq);
end
filename3 = '4bp_deletion_1_3.txt';
for ii = 281:420
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+5:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename3,seqname,new_seq);
end
filename4 = '4bp_deletion_1_4.txt';
for ii = 421:560
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+5:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename4,seqname,new_seq);
end
filename5 = '4bp_deletion_1_5.txt';
for ii = 561:700
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+5:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename5,seqname,new_seq);
end
filename6 = '4bp_deletion_1_6.txt';
for ii = 701:786
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+5:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename6,seqname,new_seq);
end


filename1 = '5bp_deletion_1_1.txt';
for ii = 1:140
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename1,seqname,new_seq);
end
filename2 = '5bp_deletion_1_2.txt';
for ii = 141:280
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename2,seqname,new_seq);
end
filename3 = '5bp_deletion_1_3.txt';
for ii = 281:420
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename3,seqname,new_seq);
end
filename4 = '5bp_deletion_1_4.txt';
for ii = 421:560
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename4,seqname,new_seq);
end
filename5 = '5bp_deletion_1_5.txt';
for ii = 561:700
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename5,seqname,new_seq);
end
filename6 = '5bp_deletion_1_6.txt';
for ii = 701:785
    new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+6:length(seq_entraper)));
    seqname = ['sequence' num2str(ii)];
    fastawrite(filename6,seqname,new_seq);
end



%%3 base deletion can also be applied, though it seems less impactful
% filename1 = '3bp_deletion_1_1.txt';
% for ii = 1:140
%     new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+4:length(seq_entraper)));
%     seqname = ['sequence' num2str(ii)];
%     fastawrite(filename1,seqname,new_seq);
% end
% filename2 = '3bp_deletion_1_2.txt';
% for ii = 141:280
%     new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+4:length(seq_entraper)));
%     seqname = ['sequence' num2str(ii)];
%     fastawrite(filename2,seqname,new_seq);
% end
% filename3 = '3bp_deletion_1_3.txt';
% for ii = 281:420
%     new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+4:length(seq_entraper)));
%     seqname = ['sequence' num2str(ii)];
%     fastawrite(filename3,seqname,new_seq);
% end
% filename4 = '3bp_deletion_1_4.txt';
% for ii = 421:560
%     new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+4:length(seq_entraper)));
%     seqname = ['sequence' num2str(ii)];
%     fastawrite(filename4,seqname,new_seq);
% end
% filename5 = '3bp_deletion_1_5.txt';
% for ii = 561:700
%     new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+4:length(seq_entraper)));
%     seqname = ['sequence' num2str(ii)];
%     fastawrite(filename5,seqname,new_seq);
% end
% filename6 = '3bp_deletion_1_6.txt';
% for ii = 701:787
%     new_seq = cat(2,seq_entraper(1:ii),seq_entraper(ii+4:length(seq_entraper)));
%     seqname = ['sequence' num2str(ii)];
%     fastawrite(filename6,seqname,new_seq);
% end_

dG5 = xlsread('5dG.xlsx');
dG4 = xlsread('4dG.xlsx');
dG5(start_1(1):start_1(1)+16) = NaN;
dG5(start_2(1):start_2(1)+19) = NaN;
dG5(start_3(1):start_3(1)+21) = NaN;
dG4(start_1(1):start_1(1)+16) = NaN;
dG4(start_2(1):start_2(1)+19) = NaN;
dG4(start_3(1):start_3(1)+21) = NaN;
start_position5 = 1:785;
start_position4 = 1:786;
plot(start_position5,dG5);
plot(start_position4,dG4);
plot(start_position5(255:270),dG5(255:270));%267
plot(start_position5(50:75),dG5(50:75));%63
plot(start_position4(15:30),dG5(15:30));%25
plot(start_position4(300:330),dG4(300:330));%319

