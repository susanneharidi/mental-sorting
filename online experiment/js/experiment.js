////////////////////////////////////////////////////////////////////////
//            JS-CODE Generic Experiment Template                     //
//            original work: CHARLEY WU                              //
////////////////////////////////////////////////////////////////////////

//TODO List when adapting this for a new experiment (not comprehensive, some things might be missing)
//1. Update Experiment parameters below
//2. Update data collectors used for the relevant beavioral data being studied
//3. Update experimental variables in assignScenario()
//4. Load relevant data in beginExperiment(), dynamically update instruction text as needed

// EXPERIMENT PARAMETERS
var totalRounds = 10, //Number of rounds
  rounds = totalRounds, //number of REMAINING rounds
  roundTracker = 0, // current round
  totalTrials = 10, //number of trials (per round)
  clicks = -1, //variable to keep track of how many remaining actions (i.e., trials) before end of round
  trialTracker = 0; //current trial number
  tracker = new Array(0), //tracker array for storing data
  investigationIndex = 0, //current trial number
  scoretotal = [], // score for calculating bonus
  scorecurrent = 0, 
  reward = 0.00, // performance based bonus
  totalEnvNumber = 20, //Number of available environments to sample from
  envOrder = getRandomSubarray(range(0,totalEnvNumber),totalRounds), //samples a random subset of environments
  currentEnv = envOrder[0], //start with first environment
  environmentList = [];
 

// data collectors for behavioral data, customize for experiment
var tscollect = [], //time stamp of choice
  xcollect= [], //variable describing choice
  ycollect= [], //variable describing choice
  zcollect = [], //reward
  initcollect = []; //initial value

// loop through each round and create an array to store data for each round
for (var i = 0; i < totalRounds; i++) {
  scoretotal[i] = 0;
  tscollect[i] = [];
  xcollect[i] = [];
  ycollect[i] = [];
  zcollect[i] = [];
  initcollect[i] = [];
};


//Declare variables not yet assigned 
var fullurl = document.location.href, //url of incoming MTurk worker, used to extract workerid and assignmentid
  scenario, //experiment scenario
  assignmentID,
  workerID,
  scenarioId,
  environment,
  experimentData,
  roundScore,
  gender=-1, //demographic info; update as needed
  age = -1;

//Access the MySQL database and retrieve a scenario id, adding the now() date time to the start time field and marking the specific scenario as completed
function assignScenario() {
  var ajaxRequest = new XMLHttpRequest();
  try{
      // Opera 8.0+, Firefox, Safari
      ajaxRequest = new XMLHttpRequest();
  } catch (e){
      // Internet Explorer Browsers
      try{
          ajaxRequest = new ActiveXObject("Msxml2.XMLHTTP");
      } catch (e) {
          try{
              ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");
          } catch (e){
              // Something went wrong
              alert("Your browser broke!");
              return false;
          }
      }
  }
  var queryString = "?action=" + 'assignScenario';
  ajaxRequest.open("GET", "databasecall.php"+queryString, false); //The main functionality of this is defined in databasecall.php
  ajaxRequest.send(null);
  var response = ajaxRequest.responseText;
  var jsonArray = JSON.parse(response);
  //TODO: Update parsed variables returned from the database, which is used to specify the experimental treatment
  scenarioId = parseInt(jsonArray['scenarioId']);
  scenario = parseInt(jsonArray['scenario']);
  
}



//Retrieve worker and assignment id from URL header, and then assigns them a scenario
function beginExperiment() { //called when participant grants consent
  // Retrieve assignmentID, workerID, ScenarioID, and environment from URL
  assignmentID = turkGetParam('assignmentId');
  workerID = turkGetParam('workerId');
  document.getElementById('MTurkID').value = workerID; //prepopulate MTurk
  //assignScenario(); //assigns a value to scenarioId and scenario, but also other variables if required
  //TODO: load relevant data for scenario assigned to participant
  //TODO: change instruction text if necessary
  //Advance the page
  clickStart('page1', 'page2');
  setButtonHandlers(); //assign functions to each of the buttons in the experiment

}



//Comprehension check questions

function instructioncheck() {
  //check if correct answers are provided
  if (document.getElementById('q1b').checked) {
    var ch1 = 1
  }
  if (document.getElementById('q2c').checked) {
    var ch2 = 1
  }
  if (document.getElementById('q3d').checked) {
    var ch3 = 1
  }
  //are all of the correct
  var checksum = ch1 + ch2 + ch3;
  if (checksum === 3) {
    //start experiment if instructions are checked correct
    clickStart('page3', 'page4');
  } 
  else {
    //if one or more answers are wrong, raise alert box
    alert('You have answered some of the questions wrong. Please try again.');
    clickStart('page3', 'page2');
  }
}

// connects first button functionality when document loads
 
document.addEventListener("DOMContentLoaded", function (event) {
  document.getElementById("buttonBeginExperiment").addEventListener("click", function () {
    beginExperiment();
  });
});


// connects all other buttons
function setButtonHandlers() {
  document.getElementById("buttonInstructions").addEventListener(clickEventType, function () {
    clickStart('page2', 'page3'); //proceed to comprehension questions
  });
  document.getElementById("buttonInstructionsCheck").addEventListener(clickEventType, function () {
    instructioncheck(); //check if comprehension questions are correct
  });
  document.getElementById("buttonGoToPageFive").addEventListener(clickEventType, function () {
    clickStart('page4', 'page5');  //start experiment
  });
  document.getElementById("nextRoundButton").addEventListener(clickEventType, function () {
    nextRound(); //go to next round
  });
  document.getElementById("goToPage7").addEventListener(clickEventType, function () {
    clickStart('page6', 'page7'); //go to demographic info
  });
  document.getElementById("finishButton").addEventListener(clickEventType, function () {
    //senddata();
    onButtonFinishPressed(); //finish experiment and save data
  });  
};



//Experiment finished
function onButtonFinishPressed() {
   //Check that all forms were filled out
  MTurkID = document.getElementById('MTurkID').value
  //gender
  if (document.getElementById('gender1').checked) {gender = 0};
  if (document.getElementById('gender2').checked) {gender = 1};
  age=document.getElementById('age').value;
  processDescription=document.getElementById('processDescription').value;
  //Check all options have been filled out
  if (MTurkID.length > 0 && age>0 && gender >-1){
    //display reward on final page
    var rewardText = "Congratulations! You earned a <b>bonus of $" + reward + "</b>, which will be automatically assigned to your MTurk account in 1-3 days once you have completed both parts of this experiment. If this is Part 1, you will be notified about Part 2 within about 24hrs.";
    change("rewardText", rewardText);
    senddata();
    }else{
      alert("Please enter your MTurk ID, age, and gender to complete the task");
    }
}


function nextRound() {
  console.log("nextRound called");
  //proceed only if there are more rounds available
  rounds = rounds - 1; // decrease remaining trials
  roundTracker = roundTracker + 1; //increase current round number
  console.log("[Debug] You have now " +rounds+ " reamaining rounds");

  if (rounds >= 0) { //still rounds remaining
    //initialize next round
    clicks = totalTrials; //reset number of clicks
    trialCounter = 0; //reset trialCounter
    //Main experiment happens here
  }
  // If remaining rounds < 0 --> game ended
  if (rounds < 0) {
    //move to final page
    clickStart('page5', 'page6');
  }
}


function debugData() {
  console.log(tscollect);
  console.log(xcollect);
  console.log(ycollect);
  console.log(zcollect);
}

function senddata() { 
  //behavioraldata history
  experimentData = {
    'tscollect': tscollect,
    'xcollect': xcollect,
    'ycollect': ycollect,
    'zcollect': zcollect,
  };

  //Initiate AJAX request
  var ajaxRequest = new XMLHttpRequest();
  try{
      // Opera 8.0+, Firefox, Safari
      ajaxRequest = new XMLHttpRequest();
  } catch (e){
    // Internet Explorer Browsers
    try{
        ajaxRequest = new ActiveXObject("Msxml2.XMLHTTP");
    } catch (e) {
        try{
            ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");
        } catch (e){
            // Something went wrong
            alert("Your browser broke!");
            return false;
        }
    }}                                                 
  var queryString = "?action=" + 'completeScenario' + '&Experiment=grid' +  '&MTurkID=' + MTurkID + '&assignmentID=' + assignmentID + '&experimentData=' + JSON.stringify(experimentData) + '&reward=' + reward + '&age=' + age + '&gender=' + gender + '&processDescription=' +processDescription  + '&scenarioId=' + scenarioId  ;
  ajaxRequest.open("GET", "databasecall.php"+queryString, false);
  ajaxRequest.send(null);
  clickStart('page7', 'page8');

}

// Touch events support
var clickEventType = "click";
window.addEventListener('touchstart', function() {
  clickEventType = "touchstart";
});


//*************UTILITIES***************************************


//changes from one page to another
function clickStart(hide, show) {
  document.getElementById(hide).style.display = "none";
  document.getElementById(show).style.display = "block";
  window.scrollTo(0, 0);
}

//changes inner HTML of div with ID=x to y
function change(x, y) {
  document.getElementById(x).innerHTML = y;
}

//adds y to inner HTML of div with ID=x
function addToDiv(x, y) {
  document.getElementById(x).innerHTML += y;
}

//Function to randomly shuffle an array:
function shuffle(o) { //v1.0
  for (var j, x, i = o.length; i; j = Math.floor(Math.random() * i), x = o[--i], o[i] = o[j], o[j] = x);
  return o;
};

//construct an array in the range [start, count]
function range(start, count) {
      return Array.apply(0, Array(count))
        .map(function (element, index) { 
          return index + start;  
      });
    };

//Randomly sample n values from an array
function getRandomSubarray(arr, size) {
  var shuffled = arr.slice(0),
    i = arr.length,
    temp, index;
  while (i--) {
    index = Math.floor((i + 1) * Math.random());
    temp = shuffled[index];
    shuffled[index] = shuffled[i];
    shuffled[i] = temp;
  }
  return shuffled.slice(0, size);
}

//load JSON file
function loadJSON(file, callback) {
  var rawFile = new XMLHttpRequest();
  rawFile.overrideMimeType("application/json");
  rawFile.open("GET", file, true);
  rawFile.onreadystatechange = function () {
    if (rawFile.readyState === 4 && rawFile.status == "200") {
      callback(rawFile.responseText);
    }
  }
  rawFile.send(null);
}

//Create normal noise distribution
function myNorm() {
  var x1, x2, rad, c;
  do {
    x1 = 2 * Math.random() - 1;
    x2 = 2 * Math.random() - 1;
    rad = x1 * x1 + x2 * x2;
  } while (rad >= 1 || rad == 0);
  c = Math.sqrt(-2 * Math.log(rad) / rad);
  return (x1 * c);
};

//average the values in an array
function average(inputArray) {
  var total = 0
  for (var i = 0; i < inputArray.length; i++) {
    total += inputArray[i];
  }
  var avg = total / inputArray.length;
  return avg;
};

//Convert cumulative score to reward value
function rewardCum(scoreTotal) {
  var r = 0,
    r_i;
  for (var i = 0; i < scoreTotal.length; i++) {
    r_i = scoreTotal[i] / (scale[i] + 5) / 190 * 1.5;
    r = r + r_i
  }
  if (r > 1.5) {
    r = 1.5; //limit to max reward, in case of any funny business
  }
  return toFixed(r, 2);
}


//single trial reward
function performanceScore(points, scale) {
  var r = 0;
  //cumulative regret (as a percentage)
  r = points / ((scale + 5) * totalTrials);
  return toFixed(r * 100);
}

function finalPerformance(scoreArray) {
  var finalScore = 0;
  for (i = 0; i < scoreArray.length; i++) { //loop through score array
    finalScore += parseInt(performanceScore(parseInt(scoreArray[i]), parseInt(scale[i])));
  }
  return toFixed(finalScore / scoreArray.length)
}
//calculate number of stars
function starsEarned(score) {
  //console.log("score: " + score);

  percentageScore = score / 100;
  //console.log("percentageScore: " + percentageScore);

  scoreOutOfFive = percentageScore * 5;
  //console.log("scoreOutOfFive: " + scoreOutOfFive);

  fixedScoreOutOfFive =  toFixed(scoreOutOfFive, 1);
  //console.log("fixedScoreOutOfFive", fixedScoreOutOfFive)

  return parseInt(fixedScoreOutOfFive) >= 5 ? 5 : fixedScoreOutOfFive;
}

function totalStarsEarned(starArray) {
  var totalStars = 0;
  for (i = 0; i < starArray.length; i++) { //loop through score array
    totalStars += parseFloat(starArray[i]);
  }
  return toFixed(totalStars)
}


//random number generator
function randomNum(min, max) {
  return Math.floor(Math.random() * (max - min + 1) + min)
}

//Display a float to a fixed percision
function toFixed(value, precision) {
  var precision = precision || 0,
    power = Math.pow(10, precision),
    absValue = Math.abs(Math.round(value * power)),
    result = (value < 0 ? '-' : '') + String(Math.floor(absValue / power));

  if (precision > 0) {
    var fraction = String(absValue % power),
      padding = new Array(Math.max(precision - fraction.length, 0) + 1).join('0');
    result += '.' + padding + fraction;
  }
  //console.log(result);
  return result;
}

// extract URL parameters (FROM: https://s3.amazonaws.com/mturk-public/externalHIT_v1.js)
function turkGetParam(name) {
  var regexS = "[\?&]" + name + "=([^&#]*)";
  var regex = new RegExp(regexS);
  var tmpURL = fullurl;
  var results = regex.exec(tmpURL);
  if (results == null) {
    return "";
  } else {
    return results[1];
  }
}

function getAge(birthDate) {

  var dob = new Date(birthDate);
  var today = new Date();
  var age2 = today.getFullYear() - dob.getFullYear();
  var m = today.getMonth() - dob.getMonth();
  if (m < 0 || (m === 0 && today.getDate() < dob.getDate())) age2--;
  return age2;
}

function getCounter(counterName) {
  if (localStorage.getItem(counterName) === null)
    return 0;
  else
    return parseInt(localStorage.getItem(counterName));
}

function incrementCounter() {
  var counter = 0;
  if (age <= 9) {
    counter = getCounter("gridsearch-counter-1");
    localStorage.setItem("gridsearch-counter-1", counter + 1);
  } else if (age > 9 && age < 18) {
    counter = getCounter("gridsearch-counter-2");
    localStorage.setItem("gridsearch-counter-2", counter + 1);
  } else if (age > 18) {
    counter = getCounter("gridsearch-counter-3");
    localStorage.setItem("gridsearch-counter-3", counter + 1);
  }
}

function expandGrid() {
    var r = [], arg = arguments, max = arg.length-1;
    function helper(arr, i) {
        for (var j=0, l=arg[i].length; j<l; j++) {
            var a = arr.slice(0); // clone arr
            a.push(arg[i][j]);
            if (i==max)
                r.push(a);
            else
                helper(a, i+1);
        }
    }
    helper([], 0);
    return r;
}

//END