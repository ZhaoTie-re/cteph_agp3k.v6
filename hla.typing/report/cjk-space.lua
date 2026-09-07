-- Drop a space that sits between two CJK characters.
--
-- Markdown needs a delimiter around a cross-reference (`@sec-x`), and inline
-- code and emphasis are usually written with spaces around them. In English
-- that space is a word separator; between two Chinese characters it is simply
-- wrong, and 22 of them survived into the rendered page even after every
-- paragraph was written on one line.
--
-- This removes the class rather than the instances: it runs after references
-- have been resolved, so it also catches the space in front of a reference that
-- expands to Chinese text ("小节 7").

local function is_cjk(cp)
  return (cp >= 0x3000 and cp <= 0x303F)      -- CJK punctuation
      or (cp >= 0x3400 and cp <= 0x4DBF)      -- extension A
      or (cp >= 0x4E00 and cp <= 0x9FFF)      -- unified ideographs
      or (cp >= 0xF900 and cp <= 0xFAFF)      -- compatibility
      or (cp >= 0xFF00 and cp <= 0xFFEF)      -- fullwidth forms
end

-- Code, Math and RawInline carry `.text`; Emph, Strong, Span, Link and Quoted
-- carry `.content`. An inline code span holding Chinese -- `基因:位置:氨基酸` --
-- is written with spaces around it and is the case this missed at first.
local function has_text(el)
  return el.t == 'Str' or el.t == 'Code' or el.t == 'Math' or el.t == 'RawInline'
end

local function last_cp(el)
  if el == nil then return nil end
  if has_text(el) then
    local cp
    for _, c in utf8.codes(el.text) do cp = c end
    return cp
  elseif el.content and #el.content > 0 then
    return last_cp(el.content[#el.content])
  end
  return nil
end

local function first_cp(el)
  if el == nil then return nil end
  if has_text(el) then
    for _, c in utf8.codes(el.text) do return c end
    return nil
  elseif el.content and #el.content > 0 then
    return first_cp(el.content[1])
  end
  return nil
end

function Inlines(inlines)
  local out = pandoc.List()
  local i = 1
  while i <= #inlines do
    local el = inlines[i]
    if el.t == 'Space' and #out > 0 and i < #inlines then
      local a, b = last_cp(out[#out]), first_cp(inlines[i + 1])
      if a and b and is_cjk(a) and is_cjk(b) then
        i = i + 1
        goto continue
      end
    end
    out:insert(el)
    i = i + 1
    ::continue::
  end
  return out
end
