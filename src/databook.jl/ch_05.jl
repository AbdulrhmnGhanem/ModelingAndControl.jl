### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ a4f3f09c-fd5f-11ee-0b84-b304f0ecd819
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
    data_path = joinpath("..", "..", "books", "databook", "DATA")
    ENV["DATADEPS_ALWAYS_ACCEPT"] = true
end;

# ╔═╡ 87218cde-7a86-4853-bd38-3ee9f1500a94
using MLDatasets,
    StatsBase,
    ImageCore,
    ImageShow,
    LinearAlgebra,
    Plots,
    MultivariateStats,
    MLJ,
    DataFrames

# ╔═╡ dc1a9923-a2a8-4ff3-a7d0-89eeb4b4d636
md"# Chapter 5 - Clustering and Classification"

# ╔═╡ 4f10506a-51b2-4251-b188-6797104778e5
md"
## Exercise 5.1
Download the MNIST data set (both training and test sets and labels) from
[http://yann.lecun.com/exdb/mnist](http://yann.lecun.com/exdb/mnist). Perform the following analysis:

1. Do an SVD analysis of the digit images. You will need to reshape each image into a column vector, and each column of your data matrix is a different image.
2. What does the singular value spectrum look like, and how many modes are necessary for good image reconstruction? (That is, what is the rank r of the digit space?)
3. What is the interpretation of the U, Σ, and V matrices?
4. On a 3D plot, project onto three selected V modes (columns) colored by their digit label, for example, columns 2, 3, and 5.

Once you have performed the above and have your data projected into PCA space, you will build a classifier to identify individual digits in the training set

5. Pick two digits. See if you can build a linear classifier (LDA) that can reasonably identify them.

6. Pick three digits. Try to build a linear classifier to identify these three now.

7. Which two digits in the data set appear to be the most difficult to separate? Quantify the accuracy of the separation with LDA on the test data.

8. Which two digits in the data set are most easy to separate? Quantify the accuracy of the separation with LDA on the test data.

9. SVM (support vector machines) and decision tree classifiers were the state of the art until about 2014. How well do these separate between all 10 digits?

10. Compare the performance between LDA, SVM, and decision trees on the hardest and easiest pair of digits to separate (from above).

Make sure to discuss the performance of your classifier on both the training and test sets.
"

# ╔═╡ b15f8541-faf1-46cc-b045-fbfdc75182a8
"""
Gavish-Donoho algorithm for unknown noise.

Adapted from [PyDMD](https://github.com/PyDMD/PyDMD/blob/d67a115805274b7862cb57d1ea4e4cd721412997/pydmd/utils.py#L19)."""
function svht(Σ::Vector, dim::Tuple{Int,Int})
    β = min(dim[1], dim[2]) / max(dim[1], dim[2])
    ω = 0.56 * β^3 - 0.95 * β^2 + 1.82 * β + 1.43
    τ = median(Σ) * ω
    rank = sum(Σ .> τ)

    if rank == 0
        @warn "SVD optimal rank is 0. The largest singular values are indistinguishable from noise. Setting rank truncation to 1."
        rank = 1
    end

    return rank
end

# ╔═╡ 1341cf2c-532c-4e06-a946-2a1fc02a7568
# ╠═╡ show_logs = false
SVMClassifier = @load SVMClassifier pkg = MLJScikitLearnInterface verbosity = 0

# ╔═╡ 9d6e2b5d-f48e-4808-ada3-e399cf6e097d
DecisionTreeClassifier = @load DecisionTreeClassifier pkg = DecisionTree verbosity = 0

# ╔═╡ d1fdbba8-ba45-488c-bf9c-766cf6dd98ae
function solve_one()
    training_data = MNIST(:train)
    X_train = reshape(training_data.features, (28 * 28, 60000))
    targets = training_data.targets

    test_data = MNIST(:test)
    X_test = reshape(test_data.features, (28 * 28, 10000))
    test_targets = test_data.targets

    # svd analysis
    U, S, Vt = svd(X_train)
    r = svht(S, size(X_train))
    V_selected = Vt'[:, [2, 3, 5]]
    projections = X_train' * V_selected

    # LDA for two and three numbers
    zs = X_train[:, findall(x -> x == 0, targets)]
    os = X_train[:, findall(x -> x == 1, targets)]
    ts = X_train[:, findall(x -> x == 2, targets)]

    targets_lda1 = vcat(zeros(size(zs, 2)), ones(size(os, 2)))
    lda1 = fit(MulticlassLDA, hcat(zs, os), targets_lda1)
    y_lda1 = MultivariateStats.predict(lda1, hcat(zs, os))

    targets_lda2 = vcat(zeros(size(zs, 2)), ones(size(os, 2)), 2 * ones(size(ts, 2)))
    lda2 = fit(MulticlassLDA, hcat(zs, os, ts), targets_lda2)
    y_lda2 = MultivariateStats.predict(lda2, hcat(zs, os, ts))

    lda3 = fit(MulticlassLDA, X_train, targets)

    X_train_lda = MultivariateStats.transform(lda3, X_train)
    X_test_lda = MultivariateStats.transform(lda3, X_test)

    classes_present = unique(targets)
    centroids = [mean(X_train_lda[:, targets.==c], dims = 2) for c in classes_present]

    function find_closest_centroid(point, centroids)
        dists = [norm(point - centroid) for centroid in centroids]
        classes_present[argmin(dists)]
    end

    y_lda3 =
        [find_closest_centroid(X_test_lda[:, i], centroids) for i = 1:size(X_test_lda, 2)]
    lda_accuracy = sum(y_lda3 .== test_targets) / length(test_targets)

    num_classes = length(centroids)
    distances = zeros(num_classes, num_classes)


    for i = 1:num_classes
        for j = i+1:num_classes
            distances[i, j] = norm(centroids[i] - centroids[j])
            distances[j, i] = distances[i, j]
        end
    end
    hard, easy = findmin(x -> x == 0 ? Inf : x, distances), findmax(distances)

    train_df = DataFrame(X_train', :auto)
    test_df = DataFrame(X_test', :auto)


    SVM = SVMClassifier(cache_size = 1000)
    svm = machine(SVM, train_df, categorical(targets))
    fit!(svm)

    svm_accuracy =
        mean(categorical(test_targets) .== coerce(MLJ.predict(svm, test_df), Multiclass))

    DecisionTree = DecisionTreeClassifier()
    decision_tree = machine(DecisionTree, train_df, categorical(targets))
    fit!(decision_tree)

    decision_tree_accuracy =
        mean(categorical(test_targets) .== mode.(MLJ.predict(decision_tree, test_df)))

    test_five = X_test[:, findall(x -> x == 5, test_targets)]
    test_three = X_test[:, findall(x -> x == 3, test_targets)]
    test_four = X_test[:, findall(x -> x == 4, test_targets)]
    test_two = X_test[:, findall(x -> x == 2, test_targets)]
    X_hard = hcat(test_three, test_five)
    hard_df = DataFrame(X_hard', :auto)
    hard_targets = vcat(5 * ones(size(test_five, 2)), 3 * ones(size(test_three, 2)))

    X_easy = hcat(test_four, test_two)
    easy_df = DataFrame(X_easy', :auto)
    easy_targets = vcat(4 * ones(size(test_four, 2)), 2 * ones(size(test_two, 2)))

    a = findall(x -> x == 5, y_lda3)
    b = findall(x -> x == 3, y_lda3)
    hard_lda_accurcy =
        mean(vcat(5 * ones(length(a)), 3 * ones(length(b))) .== test_targets[vcat(a, b)])


    a = findall(x -> x == 4, y_lda3)
    b = findall(x -> x == 2, y_lda3)
    easy_lda_accurcy =
        mean(vcat(4 * ones(length(a)), 2 * ones(length(b))) .== test_targets[vcat(a, b)])

    hard_svm_accuracy =
        mean(categorical(hard_targets) .== coerce(MLJ.predict(svm, hard_df), Multiclass))
    easy_svm_accuracy =
        mean(categorical(easy_targets) .== coerce(MLJ.predict(svm, easy_df), Multiclass))

    hard_decision_tree_accuracy =
        mean(categorical(hard_targets) .== mode.(MLJ.predict(decision_tree, hard_df)))
    easy_decision_tree_accuracy =
        mean(categorical(easy_targets) .== mode.(MLJ.predict(decision_tree, easy_df)))

    U,
    S,
    Vt,
    r,
    projections,
    targets,
    y_lda1,
    targets_lda1,
    y_lda2,
    targets_lda2,
    y_lda3,
    distances,
    easy,
    hard,
    lda_accuracy,
    svm_accuracy,
    decision_tree_accuracy,
    hard_decision_tree_accuracy,
    easy_decision_tree_accuracy,
    hard_lda_accurcy,
    easy_lda_accurcy,
    hard_svm_accuracy,
    easy_svm_accuracy
end

# ╔═╡ 2889db74-209d-4d6b-8b61-474c3e796e2c
function plot_one(sol_one)
    U,
    S,
    Vt,
    r,
    projections,
    targets,
    y_lda1,
    targets_lda1,
    y_lda2,
    targets_lda2,
    y_lda3,
    distances,
    easy,
    hard,
    lda_accuracy,
    svm_accuracy,
    decision_tree_accuracy = sol_one
    p1 = plot(S; yaxis = :log, title = "Singular values", label = missing, grid = false)
    vline!([r]; label = "optimal r")

    p2 = scatter(
        projections[:, 1],
        projections[:, 2],
        projections[:, 3],
        marker_z = targets,
        title = "3D Projection of MNIST Data",
        xlabel = "Mode 2",
        ylabel = "Mode 3",
        zlabel = "Mode 5",
    )

    p3 = plot(
        grid = false,
        axis = false,
        title = "LDA between 0 and 1",
        legend = :topleft,
        ylim = (-1, 10),
    )

    for label in [0, 1]
        points = y_lda1[:, targets_lda1.==label]'
        scatter!(p3, points, 5 * ones(length(points)), alpha = 0.5, label = label)
    end

    p4 = plot(
        grid = false,
        axis = false,
        title = "LDA between 0, 1, and 2",
        legend = :topleft,
    )
    for label in [0, 1, 2]
        points = y_lda2[:, targets_lda2.==label]
        scatter!(p4, points[1, :], points[2, :], label = label, alpha = 0.5)
    end

    p5 = scatter(
        legend = false,
        grid = false,
        aspectratio = 1,
        title = "Discrimination hardness",
        xlims = (1, 10),
    )

    # easiest and hardest to discriminate
    for i = 1:10
        for j = i:10
            scatter!(p5, [i], [j], c = :blue, markersize = 100 * distances[i, j])
        end
    end

    radius = 0.3
    theta = 0:0.01:2*pi
    plot!(
        easy[2][2] .+ radius .* cos.(theta),
        easy[2][1] .+ radius .* sin.(theta),
        linecolor = :green,
    )
    plot!(
        hard[2][2] .+ radius .* cos.(theta),
        hard[2][1] .+ radius .* sin.(theta),
        linecolor = :red,
    )

    plot(p1, p2, p3, p4, p5; layout = (3, 2))
end

# ╔═╡ 51586ab5-ae10-4274-832d-7365124f56a2
begin
    sol_one =
        UU,
        S,
        Vt,
        r,
        projections,
        targets,
        y_lda1,
        targets_lda1,
        y_lda2,
        targets_lda2,
        y_lda3,
        distances,
        easy,
        hard,
        lda_accuracy,
        svm_accuracy,
        decision_tree_accuracy,
        hard_decision_tree_accuracy,
        easy_decision_tree_accuracy,
        hard_lda_accurcy,
        easy_lda_accurcy,
        hard_svm_accuracy,
        easy_svm_accuracy = solve_one()

    plot_one(sol_one)
end

# ╔═╡ a88eaff3-9390-4357-a950-ee80e87c51e7
md"""
- The number of modes required to construct the images is $r.
##

- The easiest two digits using LDA are $(easy[2][1]) and $(easy[2][2]).
- The hardest two digits using LDA are $(hard[2][1]) and $(hard[2][2]).
##
- The $\text{accuracy}_{lda}$ = $lda_accuracy.
- The $\text{accuracy}_{svm}$ = $svm_accuracy.
- The $\text{accuracy}_{DT}$ = $decision_tree_accuracy.

##
- The $\text{accuracy}_{lda}$ (hardest) = $(round(hard_lda_accurcy; digits=4)).
- The $\text{accuracy}_{svm}$ (hardest) = $(round(hard_svm_accuracy; digits=4)).
- The $\text{accuracy}_{DT}$ (hardest) = $(round(hard_decision_tree_accuracy; digits=4)).

##
- The $\text{accuracy}_{lda}$ (easiest) = $(round(easy_lda_accurcy; digits=4)).
- The $\text{accuracy}_{svm}$ (easiest) = $(round(easy_svm_accuracy; digits=4)).
- The $\text{accuracy}_{DT}$ (easiest) = $(round(easy_decision_tree_accuracy; digits=4))


!!! remark
	Although the SVM and DT performs better on the whole testing dataset, the LDA performs much better at discirminating between the hardest two numbers: 5, and 3.
"""

# ╔═╡ 8b934310-a31c-4e13-a7f4-b74af7fbea1f
md"## Exercise 2

Download the two data sets (ORIGINAL IMAGE and CROPPED IMAGES)
from Yale Faces B. Your job is to perform an analysis of these data sets. Start with the cropped images and perform the following analysis: 

1. Do an SVD analysis of the images (where each image is reshaped into a column vector and each column is a new image).

2. What is the interpretation of the U, Σ, and V matrices?

3. What does the singular value spectrum look like and how many modes are necessary for good image reconstructions? (That is, what is the rank r of the face space?)

4. Compare the difference between the cropped (and aligned) versus uncropped images.

Face identification: see if you can build a classifier to identify individuals in the training set.

5. (Test 1) face classification: Consider the various faces and see if you can build a classifier that can reasonably identify an individual face.

6. (Test 2) gender classification: Can you build an algorithm capable of recognizing men from women?

7. (Test 3) unsupervised algorithms: In an unsupervised way, can you develop algo- rithms that automatically find patterns in the faces that naturally cluster?

(Note: You can use any (and hopefully all) of the different clustering and classification methods discussed. Be sure to compare them against each other in these tasks.)
"

# ╔═╡ 43440805-55af-4d80-ae42-a2bfc0fd0af7


# ╔═╡ Cell order:
# ╠═a4f3f09c-fd5f-11ee-0b84-b304f0ecd819
# ╠═87218cde-7a86-4853-bd38-3ee9f1500a94
# ╟─dc1a9923-a2a8-4ff3-a7d0-89eeb4b4d636
# ╟─4f10506a-51b2-4251-b188-6797104778e5
# ╠═b15f8541-faf1-46cc-b045-fbfdc75182a8
# ╠═1341cf2c-532c-4e06-a946-2a1fc02a7568
# ╠═9d6e2b5d-f48e-4808-ada3-e399cf6e097d
# ╠═d1fdbba8-ba45-488c-bf9c-766cf6dd98ae
# ╟─a88eaff3-9390-4357-a950-ee80e87c51e7
# ╟─2889db74-209d-4d6b-8b61-474c3e796e2c
# ╟─51586ab5-ae10-4274-832d-7365124f56a2
# ╟─8b934310-a31c-4e13-a7f4-b74af7fbea1f
# ╠═43440805-55af-4d80-ae42-a2bfc0fd0af7
