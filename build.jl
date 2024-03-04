using PlutoSliderServer


Export_output_dir = "./.build"
Export_cache_dir = "./.cache"
PlutoSliderServer.export_directory("./src"; Export_output_dir, Export_cache_dir)

struct Book
    title::AbstractString
    path::AbstractString
end

books = [
    Book("Practical MATLAB Modeling With Simulink", "PracticalModelingAndSimulation.jl"),
    Book("Mastering System Identification in 100 Exercises", "SysID.jl"),
    Book("Data-Driven Science and Engineering", "databook.jl"),
]

function chapters(book::Book)
    filter!(
        name -> name != "index.html",
        sort!(readdir(joinpath(Export_output_dir, book.path))),
    )
end

function format_chapter_title(title::AbstractString)
    t = replace(title, ".html" => "")
    if startswith(title, "ch_")
        t = replace(t, "ch_" => "Chapter ")
    elseif t in ["freqz"]
        t = t
    else
        t = "Exercise " * t
    end
    replace(t, "_" => ".")
end


function index_page(books::Array{Book})
    books_template = """
    <!DOCTYPE html>
    <html lang="en">
      <head>
        <meta charset="utf-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1" />

        <style>
          body {
            font-family: sans-serif;
          }

          ul > li {
            font-size: 25px;
            line-height: 50px;
          }

          a {
            text-decoration: none;
          }
        </style>

        <link rel="stylesheet" href="index.css" />
        <script src="index.js" type="module" defer></script>
      </head>
      <body>
    """

    html = books_template
    for book in books
        book_template = "<h1>$(book.title)</h1>\n<ul>\n"

        for ch in chapters(book)
            book_template *= "   <li><a href='$(book.path)/$(ch)'>$(format_chapter_title(ch))</a></li>\n"
        end
        book_template *= "</ul>\n"
        html *= book_template
        open(joinpath(Export_output_dir, book.path, "index.html"), "w") do f
            write(f, books_template * book_template * "</body>\n</html>")
        end
    end

    html *= "</body>\n</html>"

    open(joinpath(Export_output_dir, "index.html"), "w") do f
        write(f, html)
    end

end


index_page(books)
