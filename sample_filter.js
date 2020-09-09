var content = "AN02\t75\t10\n\
AN03\t351\t10\n\
AN04\t320\t10\n\
CS05\t74\t10\n\
CS07\t273\t10\n\
CS11\t312\t10\n\
CS12\t328\t10\n\
DK04\t266\t10\n\
DK05\t288\t10\n\
DK06\t67\t10\n\
DK18\t295\t10\n\
FR13\t303\t10\n\
FR19\t322\t10\n\
FR20\t91\t10\n\
FR21\t312\t10\n\
FY11\t357\t10\n\
FY13\t313\t10\n\
FY15\t377\t10\n\
HH10\t335\t10\n\
HH12\t268\t10\n\
HH16\t265\t10\n\
HH17\t316\t10\n\
JB03\t298\t10\n\
JB04\t335\t10\n\
JB05\t331\t10\n\
NB01\t248\t10\n\
NB02\t287\t10\n\
NB04\t251\t10\n\
NB16\t64\t10\n\
PA01\t286\t10\n\
PA07\t366\t10\n\
PA21\t292\t10\n\
RB05\t275\t10\n\
RB11\t305\t10\n\
RB13\t66\t10\n\
RB20\t336\t10\n\
RL08\t74\t10\n\
RL11\t229\t10\n\
RL18\t255\t10\n\
RL22\t308\t10\n\
TR02\t73\t10\n\
TR09\t286\t10\n\
TR10\t297\t10\n\
TR11\t262\t10";


console = {
    log: print,
    warn: print,
    error: print
};

var r = [];
var lines = content.split("\n");
for (var i = 0; i < lines.length; i++){
  var line = lines[i].split("\t");
  r.push({
    sampleName: line[0],
    minDP: line[2],
    maxDP: line[1]
  });
}

function record() {
  SAMPLES.forEach(function(sample) {
    r.forEach(function(sampleDat){
      if (sample==sampleDat.sampleName){
        if (sample.DP < parseInt(sampleDat.minDP)){
          sample.GT = "./."
        }
        if (sample.DP > parseInt(sampleDat.maxDP)){
          sample.GT = "./."
        }
        //console.log(sample.DP + " " + sampleDat.minDP + " " + sampleDat.maxDP + " " + sample.GT);
      }
    });
    });
}
