-- Pandoc Lua filter for the PDF build of docs/26-08-10-project-report.md.
--
-- LaTeX cannot line-break \texttt{...} content (tt hyphenation is disabled),
-- so long code-span tokens (dataset names, file paths) overflow their table
-- cells and overlap neighbouring columns. This filter re-emits every inline
-- code span as raw LaTeX with \allowbreak after the separators - _ / \ .
-- so lines may wrap at those points. Breaks are only permitted, never forced,
-- so short spans render unchanged. Code blocks are not touched.

if FORMAT ~= "latex" then
  return {}
end

function Code(el)
  local t = el.text
  -- LaTeX-escape (backslash first so inserted commands are not re-escaped)
  t = t:gsub("\\", "\\textbackslash ")
  t = t:gsub("([#$%%&_{}])", "\\%1")
  t = t:gsub("%^", "\\textasciicircum ")
  t = t:gsub("~", "\\textasciitilde ")
  -- permit (not force) a line break after each path/identifier separator
  t = t:gsub("%-", "-\\allowbreak ")
  t = t:gsub("/", "/\\allowbreak ")
  t = t:gsub("%.", ".\\allowbreak ")
  t = t:gsub("\\_", "\\_\\allowbreak ")
  t = t:gsub("\\textbackslash ", "\\textbackslash \\allowbreak ")
  return pandoc.RawInline("latex", "\\texttt{" .. t .. "}")
end
