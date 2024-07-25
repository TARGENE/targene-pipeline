
const TRAITS = [
  "J40-J47 Chronic lower respiratory diseases",
  "G25 Other extrapyramidal and movement disorders",
  "G20 Parkinson's disease",
  "L72 Follicular cysts of skin and subcutaneous tissue",
  "J34 Other disorders of nose and nasal sinuses",
  "H15-H22 Disorders of sclera, cornea, iris and ciliary body",
  "O26 Maternal care for other conditions predominantly related to pregnancy",
  "O20 Haemorrhage in early pregnancy",
  "O94-O99 Other obstetric conditions, not elsewhere classified",
  "O20-O29 Other maternal disorders predominantly related to pregnancy",
  "oesophageal disorder",
  "bacterial infection",
  "psychological/psychiatric problem",
  "heart valve problem/heart murmur",
  "other fractures",
  "glaucoma",
  "gastric/stomach ulcers",
  "mumps / epidemic parotitis",
  "other renal/kidney problem",
  "hypothyroidism/myxoedema",
  "C43 Malignant melanoma of skin",
  "D37-D48 Neoplasms of uncertain or unknown behaviour",
  "D41 Neoplasm of uncertain or unknown behaviour of urinary organs",
  "C34 Malignant neoplasm of bronchus and lung",
  "D06 Carcinoma in situ of cervix uteri",
  "D05 Carcinoma in situ of breast",
  "D04 Carcinoma in situ of skin",
  "D03 Melanoma in situ",
  "C50-C50 Malignant neoplasm of breast",
  "C81-C96 Malignant neoplasms, stated or presumed to be primary, of lymphoid, haematopoietic and related tissue",
  "Cheese intake",
  "Part of a multiple birth",
  "Ease of skin tanning",
  "Variation in diet",
  "Number of vehicles in household",
  "Skin colour",
  "Pork intake"
]

const PCS = [
  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5",
  "PC6"
]

"""
    check_fails_are_extremely_rare_traits(output, dataset)

The failures are due to traits with only very few cases (<=3) which is due to the fact that 
the dataset is quite degenerate.
"""
function check_fails_are_extremely_rare_traits(results, dataset; ncases=3)
    for result ∈ results
        Ψ = result.estimand
        @test eltype(dataset[!, Ψ.outcome]) == Bool
        nomissing = dropmissing(dataset, Symbol.([Ψ.outcome, keys(Ψ.treatment_values)...]))
        @test sum(nomissing[!, Ψ.outcome]) <= ncases
    end
end

function retrieve_failed_results(results; expected_keys=(:TMLE, :OSE))
    failed_results = (TMLE = [], OSE = [])
    for result ∈ results
        @test keys(result) == expected_keys
        @test result.TMLE isa Union{TMLE.TMLEstimate, TargetedEstimation.FailedEstimate}
        @test result.OSE isa Union{TMLE.OSEstimate, TargetedEstimation.FailedEstimate}
        if result.TMLE isa TargetedEstimation.FailedEstimate
            push!(failed_results.TMLE, result.TMLE)
        end
        if result.OSE isa TargetedEstimation.FailedEstimate
            push!(failed_results.OSE, result.OSE)
        end
    end
    return failed_results
end

function write_custom_configuration()
    config = Configuration(estimands=[
        ATE(
            outcome = "ALL",
            treatment_values = NamedTuple{(Symbol("1:238411180:T:C"), Symbol("3:3502414:T:C"))}([(control = "TT", case = "TC"), (control = "CT", case = "TT")]),
            treatment_confounders = []
        ),
      IATE(
        outcome = "ALL",
        treatment_values = NamedTuple{(Symbol("1:238411180:T:C"), Symbol("3:3502414:T:C"))}([(control = "TT", case = "TC"), (control = "CT", case = "TT")]),
        treatment_confounders = []
      ),
      CM(
        outcome = "ALL",
        treatment_values = NamedTuple{(Symbol("1:238411180:T:C"), Symbol("3:3502414:T:C"))}(["TC", "CT"]),
        treatment_confounders = []
      ),
        ATE(
        outcome = "ALL",
        treatment_values = NamedTuple{(Symbol("2:14983:G:A"),)}([(control = "AG", case = "GG")]),
        treatment_confounders = []
        ),
      CM(
        outcome = "ALL",
        treatment_values = NamedTuple{(Symbol("2:14983:G:A"),)}(["AG"]),
        treatment_confounders = []
      )
    ])
    TMLE.write_yaml(joinpath("test", "assets", "estimands.yaml"), config)
end