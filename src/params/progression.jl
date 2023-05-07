struct ProgressionParams
  # the distribution are subject to change
  dist_incubation_time::LogNormal{Float64}
  dist_symptom_onset_time::Gamma{Float64}
  dist_hospitalization_time::Gamma{Float64}
  dist_mild_recovery_time::Uniform{Float64}
  dist_death_time::LogNormal{Float64}
  dist_severity_by_age::Matrix{AliasSampler}
  hospitalization_time_ratio::Float64
end

function make_progression_params(hospitalization_time_ratio::Float64, hospitalization_multiplier::Float64, death_multiplier::Float64)
  #this data needs updating
  dist_incubation_time = LogNormal(1.23028082387, 0.47862064526)
  dist_symptom_onset_time = Gamma(0.8738003969079596, 2.9148873266517685)
  dist_hospitalization_time =  Gamma(1.1765988120148885, 2.6664347368236787)
  dist_mild_recovery_time = Uniform(11, 17)
  dist_death_time = LogNormal(1.6968381137317683, 1.2051253249941534)
  dist_severity_by_age = make_dist_severity_by_age(hospitalization_men_probs, hospitalization_women_probs, hospitalization_multiplier, death_multiplier)
  ProgressionParams(
    dist_incubation_time,
    dist_symptom_onset_time,
    dist_hospitalization_time,
    dist_mild_recovery_time,
    dist_death_time,
    dist_severity_by_age,
    hospitalization_time_ratio
  )
end

const hospitalization_women_probs = [0.00850303466539077, 0.004322657595242569, 0.0041456745792013, 0.00430996815565578, 0.0022413413513840802, 0.0021767096400654396, 0.00229799757743644, 0.00217075575105218, 0.00212196943521338, 0.00208920143602125, 0.00214577123644169, 0.0022279408083901, 0.00205643343682911, 0.00220234103686213, 0.00231011038033851, 0.00484431461593381, 0.00510275663297989, 0.00754465905562251, 0.007059174543909641, 0.00694891844548777, 0.0104332537465592, 0.0122935960831739, 0.0116229558650346, 0.0131830905366168, 0.0106756327171242, 0.007970152524288109, 0.0079062501466028, 0.0112003277266087, 0.00991144959179198, 0.0145567187410682, 0.0156873211990923, 0.0151915127217247, 0.0145286948312932, 0.0137365800582339, 0.012899438637752, 0.0124667854216415, 0.0121438432772882, 0.0161988093081708, 0.0193754664973993, 0.0193329442254486, 0.0211589692941797, 0.0215575252036099, 0.0221694951584818, 0.0278473957659549, 0.0287036295644993, 0.0317969726941198, 0.034851466943652, 0.03689307470775769, 0.0384657776904966, 0.0304449749965607, 0.0303826443862345, 0.0363034481874585, 0.0383796685057446, 0.0424155277005472, 0.0538213120805957, 0.0571178785184715, 0.0591907623661638, 0.0577040035670804, 0.0554409453002189, 0.0559547450858621, 0.0597220129173766, 0.0640461011725012, 0.0691642589894755, 0.0707368668560965, 0.0683631297561859, 0.0711670976475941, 0.0750213550577508, 0.0745655241923219, 0.0812666900509963, 0.0636436484726863, 0.0670334059317343, 0.0664673030875234, 0.0749253414774424, 0.0831350676740028, 0.0963637299998828, 0.11800000000000001, 0.135, 0.14, 0.141, 0.149336210006637, 0.152730587460395, 0.172878733473387, 0.175930847099889, 0.178453486437183, 0.204638916631492, 0.2308243468258, 0.25, 0.27, 0.29, 0.36]
const hospitalization_men_probs = [0.0080530751963399, 0.00680590446402575, 0.00612128685521843, 0.00598244068817986, 0.00604982069455727, 0.00427451969495346, 0.0043998753573024, 0.00431566917173469, 0.00204697197939929, 0.00197829367678251, 0.00197507467344147, 0.00203832046750381, 0.00212322523610964, 0.0019724466342770698, 0.0021002655778901897, 0.0022102983097089103, 0.0023, 0.00240470174797173, 0.00701141667728816, 0.0067804371254193705, 0.0067840923475894, 0.00776351002082564, 0.007428503480714841, 0.0153298481020201, 0.014629241253728, 0.0159892639896287, 0.0132468329523983, 0.0127995232961537, 0.0167420109723946, 0.0167671274651843, 0.017316121506403, 0.0199260037139624, 0.0190966406198151, 0.0183245690747419, 0.0176943736333129, 0.0162339921266343, 0.0239109942224891, 0.0256107683638328, 0.0280349561241959, 0.0303514009364933, 0.0298482734833276, 0.03936167968429329, 0.0407763114321909, 0.0409187910561636, 0.0433204120956373, 0.0447755612145209, 0.0623579428487239, 0.0485942878321935, 0.0509970702105744, 0.0538756485625393, 0.0483219725253555, 0.050923586463443, 0.05367623541625599, 0.0603094967543978, 0.0604308462478275, 0.0674672184462009, 0.0749363293267055, 0.0808881537910846, 0.0848904565952043, 0.0811149128584703, 0.0934399120110617, 0.102562265203383, 0.108125858240258, 0.124343635542474, 0.145782197166004, 0.127502444830581, 0.157372403374308, 0.160370053382243, 0.163869547940395, 0.173920086405869, 0.118509512803834, 0.11535454595601, 0.119388703835725, 0.1331419113211, 0.151680582673195, 0.169146388505085, 0.173149884792627, 0.177583462646863, 0.20495582764474302, 0.20765036752956198, 0.241720145354151, 0.284001542287219, 0.296255050053635, 0.313025655328029, 0.338762132092133, 0.35082382791373806, 0.362956157283026, 0.371758938298504, 0.39181929932390896, 0.41130052356020896, 0.481746429739453]
const max_age_hosp = length(hospitalization_men_probs)
const hospitalization_time_probs_old = [0.028939918213274615, 0.016986473733878578, 0.023592324630386914, 0.023592324630386914, 0.029569046870084933, 0.03240012582573136, 0.04687008493236867, 0.05001572821642026, 0.08776344762503932, 0.10537905001572821, 0.06511481597986789, 0.06700220195029884, 0.0978295061340044, 0.08021390374331551, 0.05536332179930796, 0.03208556149732621, 0.01918842403271469, 0.02201950298836112, 0.016042780748663103, 0.014469959106637308, 0.01384083044982699, 0.008807801195344448, 0.00817867253853413, 0.004718464926077383, 0.0053475935828877, 0.004403900597672224, 0.005662157911292859, 0.0062912865681031774, 0.0031456432840515887, 0.0031456432840515887, 0.0012582573136206354, 0.0009436929852154766, 0.0018873859704309531, 0.0012582573136206354, 0.0025165146272412707, 0.0025165146272412707, 0.0025165146272412707, 0.0015728216420257944, 0.0009436929852154766, 0.0012582573136206354, 0.0006291286568103177, 0.0012582573136206354, 0.0006291286568103177, 0.0009436929852154766, 0.0006291286568103177, 0.00031456432840515884, 0.00031456432840515884, 0.00031456432840515884, 0.00031456432840515884]
const hospitalization_time_probs = [0.09987515605493133, 0.07307532251352476, 0.08223054515189347, 0.07207657095297544, 0.06891385767790262, 0.05726175613816063, 0.05093632958801498, 0.05160216396171453, 0.04078235538909696, 0.06034124011652101, 0.1204327923429047, 0.06999583853516438, 0.04161464835622139, 0.033458177278401995, 0.02513524760715772, 0.016229712858926344, 0.01106949646275489, 0.00957136912193092, 0.00749063670411985, 0.004577611319184353, 0.0033291718684977114]
const critical_probs = [0.002472906902653, 0.0037854850012850297, 0.00402994799484427, 0.0040777513607663, 0.00675003823080385, 0.00697103847727567, 0.00722443253064095, 0.0118181678122287, 0.0127940055233372, 0.0136357361367939, 0.0140456726127444, 0.0143739319033564, 0.0165891538540999, 0.0166530394645805, 0.0170697948403793, 0.0115807911660513, 0.0119836414217112, 0.0067501338402127, 0.00774907629006527, 0.00854882255239505, 0.00707811156848143, 0.00718063342495257, 0.005818616535336579, 0.00623479608555421, 0.00720983628402203, 0.0100995161688468, 0.0115015136673871, 0.0094912789405671, 0.0111145119541469, 0.0103780570577393, 0.0103992802388452, 0.0120867139385939, 0.0141763090916946, 0.0166399138296732, 0.0201502515750186, 0.0182245590638866, 0.0197497629687048, 0.0189957537197314, 0.0189999258424321, 0.0242543128697044, 0.0222892878966428, 0.0244092944162856, 0.0272081452796555, 0.0271970195346352, 0.0297438750724428, 0.0262401345546642, 0.0333589869839148, 0.0357576930255035, 0.0384968172308294, 0.0402986379226426, 0.0440991930649188, 0.0450098179193887, 0.0463804974060005, 0.0503328734150703, 0.0481632231398294, 0.0500812118549724, 0.0533772141258145, 0.0592735475270295, 0.0699602385812662, 0.0700217744008471, 0.0730478812895055, 0.0779564399254216, 0.0779879990391362, 0.0790711126697737, 0.09884444021958001, 0.096023696753062, 0.105013691556089, 0.117513751599248, 0.12390769313845301, 0.180791771173274, 0.203571639665145, 0.226132399214654, 0.22768415072735299, 0.22852234660522402, 0.22779353838941802, 0.23464517030379298, 0.246003626468148, 0.25150570910695497, 0.252571912977583, 0.253613790932695, 0.255551315971475, 0.266459301212426, 0.287356449026597, 0.30425428704722496, 0.31579333251260305, 0.328536452437153, 0.34811541173909394, 0.361725731546196, 0.37577934164577903, 0.484142101284958]
const death_probs_age = [0.6666666666666666, 0.5945945945945946, 0.7538461538461538, 0.7801724137931034, 0.8481012658227848, 0.9402985074626866]
const hospitalization_time_sampler = AliasSampler(Int, hospitalization_time_probs)
const age_hospitalization_thresholds = Int[0, 40, 50, 60, 70, 80]


function sample_severity(rng::AbstractRNG, age::Real, gender::Bool, immunity::Bool, dist_severity_by_age::Matrix{AliasSampler})
  @assert age >= 0 "age should be non-negative"
  gender_int = gender + 1 |> UInt8
  dist = dist_severity_by_age[min(age + 1, max_age_hosp), gender_int]
  severity_int = asample(dist, rng) |> UInt8
  @assert severity_int <= 4 && severity_int > 1
  severity = severity_int |> Severity
  #reduction severity for immunited subject with some probability
  if immunity && severity_int > 1
    severity_int = 2
    severity = Mild
  end
  severity
end

function sample_if_death(rng::AbstractRNG, severity::Severity, age::Real)
  @assert age >= 0 "age should be non-negative"
  severity_int = severity |> UInt8
  if severity != Critical
    return false
  end
  death_prob = death_probs_age[agegroup(age_hospitalization_thresholds, age)]
  rand(rng) < death_prob
end

function sample_progression(rng::AbstractRNG, progression_data::ProgressionParams, age::Real, gender::Bool, immunity::Bool, time_since_immunization::Real, strain::StrainKind)
  #TODO expecting updated progression generation soon
  sample_progression(
    rng,
    age,
    gender,
    immunity,
    progression_data.dist_incubation_time,
    progression_data.dist_symptom_onset_time,
    progression_data.dist_hospitalization_time,
    progression_data.dist_mild_recovery_time,
    progression_data.dist_death_time,
    progression_data.dist_severity_by_age,
    progression_data.hospitalization_time_ratio
    
  )
end


@inline function sample_progression(rng::AbstractRNG, age::Real, gender::Bool, immunity::Bool,
    dist_incubation,
    dist_symptom_onset,
    dist_hospitalization,
    dist_mild_recovery,
    dist_death_time,
    dist_severity_by_age,
    hospitalization_time_ratio
  )

  severity = sample_severity(rng, age, gender, immunity, dist_severity_by_age)

  mild_symptoms_time = missing
  severe_symptoms_time = missing
  recovery_time = missing
  death_time = missing

  incubation_time = rand(rng, dist_incubation)

  if (severity==Mild) || (severity==Severe) || (severity==Critical)
    mild_symptoms_time = incubation_time + rand(rng, dist_symptom_onset)
  end

  severe_symptoms_time = missing
  recovery_time = missing
  death_time = missing
  if (severity==Severe) || (severity==Critical)
    severe_symptoms_time = incubation_time + rand(rng, dist_hospitalization)
    if severe_symptoms_time <= mild_symptoms_time
      mild_symptoms_time = missing
    end
    recovery_time = severe_symptoms_time + hospitalization_time_ratio * (asample(hospitalization_time_sampler, rng)-1)
  else
    if (severity==Mild)
      recovery_time = mild_symptoms_time + rand(rng, dist_mild_recovery)
    else
      if immunity # && (severity==Asymptomatic)
        recovery_time = incubation_time  # TODO if we want this # rand(rng, dist_mild_recovery)
      else # now only asymptomatic, but not vaccinated
        recovery_time = rand(rng, dist_mild_recovery)
      end
    end
  end

  if sample_if_death(rng, severity, age)
    death_time = incubation_time + rand(rng, dist_death_time)
    if death_time < severe_symptoms_time
      death_time = severe_symptoms_time 
    end
    recovery_time = missing
  end

  Progression(
    severity,
    incubation_time |> TimeDiff,
    ismissing(mild_symptoms_time) ? missing : TimeDiff(mild_symptoms_time),
    ismissing(severe_symptoms_time) ? missing : TimeDiff(severe_symptoms_time),
    ismissing(recovery_time) ? missing : TimeDiff(recovery_time),
    ismissing(death_time) ? missing : TimeDiff(death_time)
  )
end

function resample!(
  rng::AbstractRNG,
  progressions, ages::AbstractVector{T} where T <: Real,
  genders::AbstractVector{T} where T <: Bool,
  dist_incubation_time,
  dist_symptom_onset_time,
  dist_hospitalization_time,
  dist_mild_recovery_time,
  dist_death_time,
  dist_severity_by_age)

  for i in 1:length(ages)
    progressions[i] = sample_progression(rng, ages[i], genders[i], true,
      dist_incubation_time,
      dist_symptom_onset_time,
      dist_hospitalization_time,
      dist_mild_recovery_time,
      dist_death_time,
      dist_severity_by_age,
      1.0)
  end
  progressions
end

storagefloat(t::Union{Missing,TimeDiff})::Float32 = ismissing(t) ? NaN32 : Float32(t)
storagefloat(p::Progression, symbol::Symbol) = getproperty(p, symbol) |> storagefloat

function saveparams(dict, progressions::AbstractVector{Progression}, prefix::AbstractString="")
    dict[prefix*"severities"] = getproperty.(progressions, :severity)
    dict[prefix*"incubation_times"] = storagefloat.(progressions, :incubation_time)
    dict[prefix*"mild_symptoms_times"] = storagefloat.(progressions, :mild_symptoms_time)
    dict[prefix*"severe_symptoms_times"] = storagefloat.(progressions, :severe_symptoms_time)
    dict[prefix*"recovery_times"] = storagefloat.(progressions, :recovery_time)
    dict[prefix*"death_times"] = storagefloat.(progressions, :death_time)
end

function make_dist_severity_by_age(hospitalization_men_probs::Vector{T}, hospitalization_women_probs::Vector{T}, hospitalization_multiplier::Float64, death_multiplier::Float64) where T<:Real
  n = length(hospitalization_women_probs)
  hosp_man = hospitalization_men_probs * hospitalization_multiplier
  hosp_woman = hospitalization_women_probs * hospitalization_multiplier
  men_sampler = AliasSampler[]
  women_sampler = AliasSampler[]
  for idx in 1:n
    group_ids = agegroup(age_hospitalization_thresholds, idx-1)
    hosp_prob = [0, 1-hosp_man[idx], hosp_man[idx]*(1-critical_probs[idx]*death_multiplier), hosp_man[idx]*critical_probs[idx]*death_multiplier]
    push!(men_sampler, AliasSampler(Int,hosp_prob))
    hosp_prob = [0, 1-hosp_woman[idx], hosp_woman[idx]*(1-critical_probs[idx]*death_multiplier), hosp_woman[idx]*critical_probs[idx]*death_multiplier]
    push!(women_sampler, AliasSampler(Int,hosp_prob))
  end
  hcat(men_sampler,women_sampler)
end