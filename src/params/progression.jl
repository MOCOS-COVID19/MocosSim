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

const hospitalization_women_probs = [0.01313763700843207, 0.008503034665390773, 0.0043226575952425746,0.004145674579201299, 0.004309968155655776, 0.0022413413513840776,0.0021767096400654414, 0.002297997577436438, 0.0021707557510521774,0.0021219694352133836, 0.004043488439064005, 0.0037925493115250637,0.0038459480084456843, 0.0020564334368291057,0.0022023410368621314, 0.00231011038033851, 0.004844314615933812,0.005102756632979894, 0.007544659055622508, 0.007059174543909641,0.006948918445487766, 0.010433253746559228, 0.012293596083173892,0.011622955865034574, 0.013183090536616813, 0.010675632717124234,0.007970152524288112, 0.007906250146602801, 0.011200327726608663,0.009911449591791981, 0.014556718741068217, 0.015687321199092316,0.015191512721724712, 0.014528694831293218, 0.013736580058233853,0.01289943863775203, 0.012466785421641547, 0.012143843277288173,0.016198809308170827, 0.019375466497399275, 0.019332944225448557,0.021158969294179736, 0.021557525203609894, 0.022169495158481774,0.027847395765954894, 0.028703629564499317, 0.031796972694119796,0.034851466943652015, 0.036893074707757684, 0.03846577769049655,0.030444974996560667, 0.030382644386234534, 0.03630344818745851,0.03837966850574455, 0.04241552770054717, 0.05382131208059565,0.05711787851847152, 0.05919076236616383, 0.057704003567080414,0.055440945300218925, 0.060139689637103394, 0.05972201291737657,0.06404610117250119, 0.06916425898947548, 0.07073686685609654,0.06836312975618591, 0.07116709764759413, 0.0750213550577508,0.0745655241923219, 0.08126669005099627, 0.06364364847268633,0.06703340593173425, 0.06646730308752338, 0.07492534147744241,0.08313506767400276, 0.09636372999988278, 0.09331193061138679,0.08374570691045126, 0.08484591325944799, 0.07976811494737908,0.149336210006637, 0.15273058746039497, 0.1728787334733867,0.1759308470998886, 0.17845348643718326, 0.23082434682579953]
const hospitalization_men_probs = [0.0080530751963399, 0.0068059044640257525, 0.006121286855218429,0.005982440688179858, 0.00604982069455727, 0.004274519694953461,0.0062286247000953405, 0.004315669171734686, 0.002046971979399289,0.0019782936767825106, 0.0038474597804615835,0.0018451560289357049, 0.0018272441222453516, 0.001972446634277074,0.002100265577890194, 0.002210298309708905, 0.0023, 0.002404701747971731, 0.007011416677288157, 0.006780437125419366,0.0067840923475894046, 0.0077635100208256385, 0.007428503480714841,0.015329848102020118, 0.01462924125372798, 0.015989263989628688,0.013246832952398258, 0.012799523296153727, 0.01674201097239455,0.016767127465184325, 0.01731612150640295, 0.019926003713962446,0.01909664061981508, 0.018324569074741875, 0.017694373633312926,0.016233992126634254, 0.02391099422248912, 0.025610768363832807,0.028034956124195917, 0.030351400936493314, 0.02984827348332758,0.039361679684293306, 0.040776311432190925, 0.04091879105616364,0.04332041209563725, 0.04477556121452089, 0.062357942848723924,0.04859428783219349, 0.05099707021057441, 0.05387564856253928,0.04832197252535554, 0.05092358646344297, 0.05912515974731158,0.06030949675439781, 0.06043084624782745, 0.07035289093645591,0.08165284926012069, 0.08088815379108455, 0.08489045659520425,0.08111491285847028, 0.09343991201106173, 0.1025622652033834,0.10812585824025751, 0.14281116167416555, 0.14578219716600413,0.12750244483058099, 0.15737240337430794, 0.16037005338224267,0.16386954794039457, 0.17392008640586887, 0.11850951280383412,0.11535454595601041, 0.11938870383572489, 0.1331419113211,0.15168058267319467, 0.16914638850508523, 0.14961181918044905,0.16699798887857514, 0.16699995004059387, 0.16369672916907313,0.25279957912113465, 0.2840015422872193, 0.2962550500536351,0.3130256553280294, 0.338762132092133, 0.3629561572830264]
const max_age_hosp = length(hospitalization_men_probs)
const hospitalization_time_probs_old = [0.028939918213274615, 0.016986473733878578, 0.023592324630386914, 0.023592324630386914, 0.029569046870084933, 0.03240012582573136, 0.04687008493236867, 0.05001572821642026, 0.08776344762503932, 0.10537905001572821, 0.06511481597986789, 0.06700220195029884, 0.0978295061340044, 0.08021390374331551, 0.05536332179930796, 0.03208556149732621, 0.01918842403271469, 0.02201950298836112, 0.016042780748663103, 0.014469959106637308, 0.01384083044982699, 0.008807801195344448, 0.00817867253853413, 0.004718464926077383, 0.0053475935828877, 0.004403900597672224, 0.005662157911292859, 0.0062912865681031774, 0.0031456432840515887, 0.0031456432840515887, 0.0012582573136206354, 0.0009436929852154766, 0.0018873859704309531, 0.0012582573136206354, 0.0025165146272412707, 0.0025165146272412707, 0.0025165146272412707, 0.0015728216420257944, 0.0009436929852154766, 0.0012582573136206354, 0.0006291286568103177, 0.0012582573136206354, 0.0006291286568103177, 0.0009436929852154766, 0.0006291286568103177, 0.00031456432840515884, 0.00031456432840515884, 0.00031456432840515884, 0.00031456432840515884]
const hospitalization_time_probs = [0.09987515605493133, 0.07307532251352476, 0.08223054515189347, 0.07207657095297544, 0.06891385767790262, 0.05726175613816063, 0.05093632958801498, 0.05160216396171453, 0.04078235538909696, 0.06034124011652101, 0.1204327923429047, 0.06999583853516438, 0.04161464835622139, 0.033458177278401995, 0.02513524760715772, 0.016229712858926344, 0.01106949646275489, 0.00957136912193092, 0.00749063670411985, 0.004577611319184353, 0.0033291718684977114]
const critical_probs = [0.030927835051546393, 0.04161979752530934, 0.050980392156862744, 0.09044834307992203, 0.12330905306971904, 0.17335058214747737]
const death_probs_age = [0.6666666666666666, 0.5945945945945946, 0.7538461538461538, 0.7801724137931034, 0.8481012658227848, 0.9402985074626866]
const hospitalization_time_sampler = AliasSampler(Int, hospitalization_time_probs)
const age_hospitalization_thresholds = Int[0, 40, 50, 60, 70, 80]


function sample_severity(rng::AbstractRNG, age::Real, gender::Bool, immunity::ImmunityState, dist_severity_by_age::Matrix{AliasSampler})
  @assert age >= 0 "age should be non-negative"
  gender_int = gender + 1 |> UInt8
  dist = dist_severity_by_age[min(age + 1, max_age_hosp), gender_int]
  severity_int = asample(dist, rng) |> UInt8
  @assert severity_int <= 4 && severity_int > 1
  severity = severity_int |> Severity
  #reduction severity for immunited subject with some probability
  if immunity == Immunity && severity_int > 1
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

function sample_progression(rng::AbstractRNG, progression_data::ProgressionParams, age::Real, gender::Bool, immunity::ImmunityState, time_since_immunization::Real, strain::StrainKind)
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


@inline function sample_progression(rng::AbstractRNG, age::Real, gender::Bool, immunity::ImmunityState,
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
      if immunited(immunity) # && (severity==Asymptomatic)
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
    progressions[i] = sample_progression(rng, ages[i], genders[i],
      dist_incubation_time,
      dist_symptom_onset_time,
      dist_hospitalization_time,
      dist_mild_recovery_time,
      dist_death_time,
      dist_severity_by_age)
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
    hosp_prob = [0, 1-hosp_man[idx], hosp_man[idx]*(1-critical_probs[group_ids]*death_multiplier), hosp_man[idx]*critical_probs[group_ids]*death_multiplier] 
    push!(men_sampler, AliasSampler(Int,hosp_prob))
    hosp_prob = [0, 1-hosp_woman[idx], hosp_woman[idx]*(1-critical_probs[group_ids]*death_multiplier), hosp_woman[idx]*critical_probs[group_ids]*death_multiplier] 
    push!(women_sampler, AliasSampler(Int,hosp_prob))
  end
  hcat(men_sampler,women_sampler)
end