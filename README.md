# [Title]	

[Authors, Year]						

## Data and Code Used for Statistical Analyses

### The Task 

- task comprised a speeded go/no-go task with an embedded word categorization task<br><br>
- go/no-go task
	- a white arrow pointing either upward or downward either turned 
	- (a) green and kept its initial orientation (two-thirds of the trials; go trials)
	- (b) turquoise and kept its orientation (one-sixth of the trials; no-go trials), or 
	- (c) green but reversed its orientation (one-sixth of the trials; no-go trials)<br><br>
- four response types in the go/no-go task were differentiated: 
	- SH: slow hits (i.e., correct responses in go trials above the individual RT limit)
	- FH: fast hits (i.e., correct responses in go trials below the RT limit)
	- FA: false alarms (i.e., erroneous responses in no-go trials)
	- IR: correctly inhibited responses (i.e., successful response inhibitions in no-go trials)<br><br>
- word categorization task
	- after each go/no-go trial, an affective word was presented, which the participants were asked to categorize as either positive or negative<br><br>
- during the task, the skin conductance response (SCR) was recorded 
	- for statistical analyses, the integrated SCR (ISCR; Benedek & Kaernbach, 2010) was extracted within a time window from 1-4 s after the participant's response in the go/no-go task<br><br>

### The Design

- 4 x 2 design: go/no-go response type (SH, FH, FA, IR) and word valence (positive, negative) as within-participants factors<br><br>

### The Data Set ("Single_Trial_Data.rda")

- 15480 observations, 16 variables (30 participants, 516 trials per participant)

| VARIABLE          	| DESCRIPTION                        	| VALUES                                                                                                                                                                 	|
|-------------------	|------------------------------------	|------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| subject_id        	| Participant identifier             	| 1â€“30                                                                                                                                                                   	|
| trial             	| Trial number within the task       	| 1â€“516 per participant                                                                                                                                                  	|
| block             	| Experimental block                 	| 1â€“9 per participant                                                                                                                                                    	|
| gng_response_type 	| Response type in the go/no-go task 	| SH = slow hit <br> FH = fast hit <br>  FA = false alarm <br> IR = inhibited response <br> miss = missing response in go trial <br> miss = missing response in go trial 	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|
|                   	|                                    	|                                                                                                                                                                        	|                                                                                                                                                                  	|
_____________________________________________________________________________________________________________________________________________________________

VARIABLE			DESCRIPTION							VALUES
_____________________________________________________________________________________________________________________________________________________________

subject_id			Participant identifier 					1-30 
_____________________________________________________________________________________________________________________________________________________________

trial				Trial number within the task					1-516 per participant
_____________________________________________________________________________________________________________________________________________________________

block				Experimental block 						1-9 per participant 
_____________________________________________________________________________________________________________________________________________________________

gng_response_type		Response type in the go/no-go task				SH = slow hit
												FH = fast hit
												FA = false alarm
												IR = inhibited response
												miss = missing response in go trial
												false_key = response made with a word categorization key
_____________________________________________________________________________________________________________________________________________________________

gng_rt				Response time (RT) in the go/no-go task				Measured in milliseconds
												NAs for trials in which no response was made (IR, miss)
_____________________________________________________________________________________________________________________________________________________________

gng_rt_invalid			Indication whether RT in the go/no-go task was < 100 ms  	TRUE 
				or > 700 ms							FALSE
_____________________________________________________________________________________________________________________________________________________________
													
gng_rt_inverse			Inverse-transformed RT in the go/no-go task			-1000/RT in milliseconds
_____________________________________________________________________________________________________________________________________________________________

word				Affective word presented in the word-categorization task	60 German words (30 positive, 30 negative) from the Berlin 
												Affective Word List Reloaded (Vo et al., 2009)
_____________________________________________________________________________________________________________________________________________________________

word_valence			Valence of the affective word					pos = positive
												neg = negative
_____________________________________________________________________________________________________________________________________________________________
													
condition			Specific trial type (combination of word_valence and 		12 conditions:
				gng_response_type)						neg_after_[SH/FH/FA/IR/miss/false_key]
												neg_after_[SH/FH/FA/IR/miss/false_key]
												(conditions containing miss and false_key of no interest)
_____________________________________________________________________________________________________________________________________________________________

word_accuracy			Correctness of the word categorization				correct = correct categorization
												incorrect = incorrect categorization
												miss = missing response
												false_key = response made with the go/no-go key
_____________________________________________________________________________________________________________________________________________________________

word_rt				Response time in the word categorization task			Measured in milliseconds
												NAs for trials in which no response was made (miss)
_____________________________________________________________________________________________________________________________________________________________

word_rt_outlier			Indication whether RT in the word categorization task 		TRUE 
				was > 3 median absolute deviations (MAD) above or below  	FALSE 
				a participant's median RT computed per condition
_____________________________________________________________________________________________________________________________________________________________

word_rt_inverse			Inverse-transformed RT in the word categorization task		-1000/RT in milliseconds
_____________________________________________________________________________________________________________________________________________________________

iscr_gng_resp			Integrated SCR extracted within a time window from 1–4 s 	 Measured in microsiemens
				after the response in the go/no-go task
_____________________________________________________________________________________________________________________________________________________________

trial_followed_or_preceded_	Indication whether trial was followed/preceded by any incor-	TRUE 
by_any_incorr_resp		rect response in the go/no-go or the word categorization task	FALSE 
_____________________________________________________________________________________________________________________________________________________________



Analyses were conducted with R version 3.6.1 (2019-07-05) using R Studio Version 1.2.5001 © 2009-2019 RStudio, Inc. on a Windows computer.



Benedek, M., & Kaernbach, C. (2010). A continuous measure of phasic electrodermal activity. Journal of Neuroscience Methods, 190(1), 80-91. 
	https://doi.org/10.1016/j.jneumeth.2010.04.028 
	
Vo, M. L., Conrad, M., Kuchinke, L., Urton, K., Hofmann, M. J., & Jacobs, A. M. (2009). The Berlin Affective Word List Reloaded (BAWL-R). Behavior Research 
	Methods, 41(2), 534-538. https://doi.org/10.3758/BRM.41.2.534 	

