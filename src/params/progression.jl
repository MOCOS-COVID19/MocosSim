struct Progression
    severity::Severity
    # times given with respect to the infection time
    incubation_time::TimeDiff
    mild_symptoms_time::Union{Missing,TimeDiff}
    severe_symptoms_time::Union{Missing,TimeDiff}
    critical_symptoms_time::Union{Missing,TimeDiff}
    recovery_time::Union{Missing,TimeDiff}
    death_time::Union{Missing, TimeDiff}
    #Progression(severity::Severity, incubation_time::Real, mild_time::Real, severe_time::Real, recovery_time) = incubation < mild_time < severe_time < recovery_time
end

Progression() = Progression(Asymptomatic, missing, missing, missing, missing, missing)

# Asymptotic not allowed as yet (some asserts prevent them)
const severity_dists = [
  Categorical([0,  0.852,  0.144,  0.004]),
  Categorical([0,  0.848,  0.144,  0.008]),
  Categorical([0,  0.832,  0.141,  0.027]),
  Categorical([0,  0.793,  0.134,  0.073]),
  Categorical([0,  0.716,  0.121,  0.163]),
  Categorical([0,  0.596,  0.102,  0.302]),
]

const hospitalization_women_probs = [0.01313763700843207, 0.008503034665390773, 0.0043226575952425746,0.004145674579201299, 0.004309968155655776, 0.0022413413513840776,0.0021767096400654414, 0.002297997577436438, 0.0021707557510521774,0.0021219694352133836, 0.004043488439064005, 0.0037925493115250637,0.0038459480084456843, 0.0020564334368291057,0.0022023410368621314, 0.00231011038033851, 0.004844314615933812,0.005102756632979894, 0.007544659055622508, 0.007059174543909641,0.006948918445487766, 0.010433253746559228, 0.012293596083173892,0.011622955865034574, 0.013183090536616813, 0.010675632717124234,0.007970152524288112, 0.007906250146602801, 0.011200327726608663,0.009911449591791981, 0.014556718741068217, 0.015687321199092316,0.015191512721724712, 0.014528694831293218, 0.013736580058233853,0.01289943863775203, 0.012466785421641547, 0.012143843277288173,0.016198809308170827, 0.019375466497399275, 0.019332944225448557,0.021158969294179736, 0.021557525203609894, 0.022169495158481774,0.027847395765954894, 0.028703629564499317, 0.031796972694119796,0.034851466943652015, 0.036893074707757684, 0.03846577769049655,0.030444974996560667, 0.030382644386234534, 0.03630344818745851,0.03837966850574455, 0.04241552770054717, 0.05382131208059565,0.05711787851847152, 0.05919076236616383, 0.057704003567080414,0.055440945300218925, 0.060139689637103394, 0.05972201291737657,0.06404610117250119, 0.06916425898947548, 0.07073686685609654,0.06836312975618591, 0.07116709764759413, 0.0750213550577508,0.0745655241923219, 0.08126669005099627, 0.06364364847268633,0.06703340593173425, 0.06646730308752338, 0.07492534147744241,0.08313506767400276, 0.09636372999988278, 0.09331193061138679,0.08374570691045126, 0.08484591325944799, 0.07976811494737908,0.149336210006637, 0.15273058746039497, 0.1728787334733867,0.1759308470998886, 0.17845348643718326, 0.23082434682579953]
const hospitalization_men_probs = [0.0080530751963399, 0.0068059044640257525, 0.006121286855218429,0.005982440688179858, 0.00604982069455727, 0.004274519694953461,0.0062286247000953405, 0.004315669171734686, 0.002046971979399289,0.0019782936767825106, 0.0038474597804615835,0.0018451560289357049, 0.0018272441222453516, 0.001972446634277074,0.002100265577890194, 0.002210298309708905, 0.0023, 0.002404701747971731, 0.007011416677288157, 0.006780437125419366,0.0067840923475894046, 0.0077635100208256385, 0.007428503480714841,0.015329848102020118, 0.01462924125372798, 0.015989263989628688,0.013246832952398258, 0.012799523296153727, 0.01674201097239455,0.016767127465184325, 0.01731612150640295, 0.019926003713962446,0.01909664061981508, 0.018324569074741875, 0.017694373633312926,0.016233992126634254, 0.02391099422248912, 0.025610768363832807,0.028034956124195917, 0.030351400936493314, 0.02984827348332758,0.039361679684293306, 0.040776311432190925, 0.04091879105616364,0.04332041209563725, 0.04477556121452089, 0.062357942848723924,0.04859428783219349, 0.05099707021057441, 0.05387564856253928,0.04832197252535554, 0.05092358646344297, 0.05912515974731158,0.06030949675439781, 0.06043084624782745, 0.07035289093645591,0.08165284926012069, 0.08088815379108455, 0.08489045659520425,0.08111491285847028, 0.09343991201106173, 0.1025622652033834,0.10812585824025751, 0.14281116167416555, 0.14578219716600413,0.12750244483058099, 0.15737240337430794, 0.16037005338224267,0.16386954794039457, 0.17392008640586887, 0.11850951280383412,0.11535454595601041, 0.11938870383572489, 0.1331419113211,0.15168058267319467, 0.16914638850508523, 0.14961181918044905,0.16699798887857514, 0.16699995004059387, 0.16369672916907313,0.25279957912113465, 0.2840015422872193, 0.2962550500536351,0.3130256553280294, 0.338762132092133, 0.3629561572830264]
const max_age_hosp = length(hospitalization_men_probs)
const hospitalization_time_probs = [0.02317291, 0.02139037, 0.02558457, 0.0285205, 0.03627975, 0.04309531, 0.06154975, 0.08105274, 0.08608577, 0.07916536, 0.07664884, 0.08168187, 0.07780224, 0.0558876 , 0.03554577, 0.02443116, 0.01908357, 0.01751075, 0.01478452, 0.01237286, 0.01027577, 0.00723498, 0.00608158, 0.00482332, 0.00513788, 0.00545245, 0.00503303, 0.00419419, 0.00251651, 0.00178253, 0.00136311, 0.00136311, 0.00188739, 0.0020971 , 0.00251651, 0.00220195, 0.00167768, 0.00125826, 0.00094369, 0.00104855, 0.00083884, 0.00094369, 0.00073398, 0.00062913, 0.00041942, 0.00031456, 0.00031456, 0.00031456, 0.00031456]

const critical_probs = 2 * [0.030927835051546393, 0.04161979752530934, 0.050980392156862744, 0.09044834307992203, 0.12330905306971904, 0.17335058214747737]
const death_probs_age = [0.6666666666666666, 0.5945945945945946, 0.7538461538461538, 0.7801724137931034, 0.8481012658227848, 0.9402985074626866]
const hospitalization_time_sampler = AliasSampler(Int, hospitalization_time_probs)
const age_hospitalization_thresholds = Int[0, 40, 50, 60, 70, 80]

const vaccination_severe_effectiveness = 0.875
const vaccination_critical_effectiveness = 0.92
const vaccination_mild_effectiveness = 0.33
#[0, 0.6, 0.85, 0.85] # Asymptomatic=1 Mild Severe Critical
function sample_severity(rng::AbstractRNG, age::Real, gender::Bool, severity_dists_ages, vaccinated::Bool)
  if age < 0
    error("age should be non-negative")
  end
  group_ids = agegroup(age_hospitalization_thresholds, age)
  dist = severity_dists[group_ids]
  idx = gender + 1 |> Int32
  #prob = hospitalization_probs[idx][age < max_age_hosp ? age + 1 : max_age_hosp]
  #dist = Categorical([0,  1.0 - prob,  prob *(1-critical_probs[group_ids]) ,  prob *critical_probs[group_ids]])
  dist = severity_dists_ages[idx][age < max_age_hosp ? age + 1 : max_age_hosp]
  severity_int = rand(rng, dist)
  severity = rand(rng, dist) |> Severity
  if vaccinated
    severity = Asymptomatic
  end
  severity
end

function sample_if_death(rng::AbstractRNG, severity::Severity, age::Real)
  if age < 0
    error("age should be non-negative")
  end
  severity_int = severity |> UInt8
  if severity != Critical
    return false
  end
  death_prob = death_probs_age[agegroup(age_hospitalization_thresholds, age)]
  rand(rng) < death_prob
end

@inline function sample_progression(rng::AbstractRNG, age::Real, gender::Bool,
    dist_incubation,
    dist_symptom_onset,
    dist_hospitalization,
    dist_mild_recovery,
    dist_severe_recovery,
    dist_death_time,
    severity_dists_ages,
    age_vaccination_thresholds,
    vaccination_uptakes_probs_age
  )
  @assert length(vaccination_uptakes_probs_age) == length(age_vaccination_thresholds)
  vaccine_prob = vaccination_uptakes_probs_age[agegroup(age_vaccination_thresholds, age)]
  vaccinated = rand(rng) < vaccine_prob

  severity = sample_severity(rng, age, gender, severity_dists_ages, vaccinated)

  mild_symptoms_time = missing
  severe_symptoms_time = missing
  recovery_time = missing
  death_time = missing

  incubation_time = rand(rng, dist_incubation)

  if (severity==Mild) || (severity==Severe) || (severity==Critical)
    mild_symptoms_time = incubation_time + rand(rng, dist_symptom_onset)
  end
  severe_symptoms_time = missing
  critical_symptoms_time = missing
  recovery_time = missing
  death_time = missing
  if (severity==Severe) || (severity==Critical)
    severe_symptoms_time = incubation_time + rand(rng, dist_hospitalization)
    if severe_symptoms_time <= mild_symptoms_time
      mild_symptoms_time = missing
    end
    if (severity==Critical)
      critical_symptoms_time = severe_symptoms_time
    end
    recovery_time = severe_symptoms_time + asample(hospitalization_time_sampler, rng)
  else
    if (severity==Mild)
      recovery_time = mild_symptoms_time + rand(rng, dist_mild_recovery)
    elseif vaccinated  # && (severity==Asymptomatic)
      recovery_time = incubation_time
    else # now only asymptomatic, but not vaccinated
      recovery_time = rand(rng, dist_mild_recovery)
    end
  end

  if sample_if_death(rng, severity, age) # vaccine effect was already considered, no need to modify here
    death_time = incubation_time + rand(rng, dist_death_time)
    if death_time < severe_symptoms_time
      severe_symptoms_time = death_time
      if isequal(mild_symptoms_time, missing)
        mild_symptoms_time = missing
      elseif severe_symptoms_time <= mild_symptoms_time
        mild_symptoms_time = missing
      end
    end
    recovery_time = missing
  end

  Progression(
    severity,
    incubation_time |> TimeDiff,
    ismissing(mild_symptoms_time) ? missing : TimeDiff(mild_symptoms_time),
    ismissing(severe_symptoms_time) ? missing : TimeDiff(severe_symptoms_time),
    ismissing(critical_symptoms_time) ? missing : TimeDiff(critical_symptoms_time),
    ismissing(recovery_time) ? missing : TimeDiff(recovery_time),
    ismissing(death_time) ? missing : TimeDiff(death_time)
  )
end

function resample!(
  rng::AbstractRNG,
  progressions, ages::AbstractVector{T} where T <: Real,
  genders::AbstractVector{Bool},
  dist_incubation_time,
  dist_symptom_onset_time,
  dist_hospitalization_time,
  dist_mild_recovery_time,
  dist_severe_recovery_time,
  dist_death_time,
  severity_dists_ages,
  age_vaccination_thresholds,
  vaccination_uptakes_probs_age)

  for i in 1:length(ages) # ages is +1 as we have older population now
    progressions[i] = sample_progression(rng, ages[i] + 1, genders[i],
      dist_incubation_time,
      dist_symptom_onset_time,
      dist_hospitalization_time,
      dist_mild_recovery_time,
      dist_severe_recovery_time,
      dist_death_time,
      severity_dists_ages,
      age_vaccination_thresholds,
      vaccination_uptakes_probs_age)
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
    dict[prefix*"critical_symptoms_times"] = storagefloat.(progressions, :critical_symptoms_time)
    dict[prefix*"recovery_times"] = storagefloat.(progressions, :recovery_time)
    dict[prefix*"death_times"] = storagefloat.(progressions, :death_time)
end

function make_severity_dists_ages(hospitalization_men_probs::Vector{T}, hospitalization_women_probs::Vector{T}) where T<:Real
  n = length(hospitalization_women_probs)
  severity_dists_ages = [Vector{Categorical}(undef, n),Vector{Categorical}(undef, n)]
  for idx in 1:n
    group_ids = agegroup(age_hospitalization_thresholds, idx-1)
    severity_dists_ages[1][idx] = Categorical([0,  1-hospitalization_men_probs[idx],  hospitalization_men_probs[idx]*(1-critical_probs[group_ids]),  hospitalization_men_probs[idx]*critical_probs[group_ids]])
    severity_dists_ages[2][idx] = Categorical([0,  1-hospitalization_women_probs[idx],  hospitalization_women_probs[idx]*(1-critical_probs[group_ids]),  hospitalization_women_probs[idx]*critical_probs[group_ids]])
  end
  severity_dists_ages
end