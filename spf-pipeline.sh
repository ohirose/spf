DIRCNF=./setting/conf
DIROUT=./output/data
DIRBIN=./bin

## If you'd like to skip some of the following steps,
## comment out the steps you'd like to skip by attaching
## '#' symbols to the heads of the corresponding sentences
${DIRBIN}/spf-convert ${DIRCNF}-convert.txt ${DIROUT} 
${DIRBIN}/spf-detect  ${DIRCNF}-detect.txt  ${DIROUT}
${DIRBIN}/spf-track   ${DIRCNF}-track.txt   ${DIROUT}
${DIRBIN}/spf-view    ${DIRCNF}-view.txt    ${DIROUT} ${DIROUT}.trk

