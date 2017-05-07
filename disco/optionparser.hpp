/*
 *  BEANDisco: command line option parsing, extends boost/program_options
 *  
 *  Copyright 2015 Teppo Niinimäki <teppo.niinimaki(at)helsinki.fi>
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <boost/program_options.hpp>


#ifndef OPTIONPARSER_HPP
#define OPTIONPARSER_HPP

template <typename T>
T validate(const std::string& s) {
	try {
		return boost::lexical_cast<unsigned int>(s);
	}
	catch (const boost::bad_lexical_cast&) {
		throw boost::program_options::invalid_option_value(s);
	}
}


namespace boost { namespace program_options {

template<class T, class charT = char>
class typed_value_extended : public boost::program_options::value_semantic_codecvt_helper<charT>,
					public boost::program_options::typed_value_base
{
public:
	using validate_function = std::function<void (boost::any& v, 
			const std::vector<std::basic_string<charT> >& s)>;

	typed_value_extended(T* store_to) 
	: m_store_to(store_to), m_composing(false),
	  m_multitoken(false), m_zero_tokens(false),
	  m_required(false)
	{} 

	typed_value_extended* default_value(const T& v)
	{
		m_default_value = boost::any(v);
		m_default_value_as_text = boost::lexical_cast<std::string>(v);
		return this;
	}

	typed_value_extended* default_value(const T& v, const std::string& textual)
	{
		m_default_value = boost::any(v);
		m_default_value_as_text = textual;
		return this;
	}

	typed_value_extended* default_value_text(const std::string& textual)
	{
		m_default_value_as_text = boost::lexical_cast<std::string>(textual);
		return this;
	}

	typed_value_extended* implicit_value(const T &v)
	{
		m_implicit_value = boost::any(v);
		m_implicit_value_as_text =
			boost::lexical_cast<std::string>(v);
		return this;
	}

	typed_value_extended* value_name(const std::string& name)
	{
		m_value_name = name;
		return this;
	}

	typed_value_extended* implicit_value(const T &v, const std::string& textual)
	{
		m_implicit_value = boost::any(v);
		m_implicit_value_as_text = textual;
		return this;
	}


	typed_value_extended* validator(validate_function f)
	{
		m_validate = f;
		return this;
	}

	typed_value_extended* validator(std::function<void (boost::any&, const std::string&)> f)
	{
		m_validate = [f] (boost::any& v, const std::vector<std::basic_string<charT> >& s) {
			validators::check_first_occurrence(v);
			const std::string& ss = validators::get_single_string(s);
			f(v, ss);
		};
		return this;
	}

	typed_value_extended* validator(std::function<T (const std::string&)> f)
	{
		m_validate = [f] (boost::any& v, const std::vector<std::basic_string<charT> >& s) {
			validators::check_first_occurrence(v);
			const std::string& ss = validators::get_single_string(s);
			v = boost::any(f(ss));
		};
		return this;
	}


	typed_value_extended* notifier(boost::function1<void, const T&> f)
	{
		m_notifier = f;
		return this;
	}

	typed_value_extended* composing()
	{
		m_composing = true;
		return this;
	}

	typed_value_extended* multitoken()
	{
		m_multitoken = true;
		return this;
	}

	typed_value_extended* zero_tokens() 
	{
		m_zero_tokens = true;
		return this;
	}
		
	typed_value_extended* required()
	{
		m_required = true;
		return this;
	}

public: // value semantic overrides

	std::string name() const
	{
		std::string const& var = (m_value_name.empty() ? arg : m_value_name);
		if (!m_implicit_value.empty() && !m_implicit_value_as_text.empty()) {
			std::string msg = "[=" + var + "(=" + m_implicit_value_as_text + ")]";
			if (!m_default_value_as_text.empty())
				msg += " (=" + m_default_value_as_text + ")";
			return msg;
		}
		else if (!m_default_value_as_text.empty()) {
			return var + " (=" + m_default_value_as_text + ")";
		} else {
			return var;
		}
	}

	bool is_composing() const { return m_composing; }

	unsigned min_tokens() const
	{
		if (m_zero_tokens || !m_implicit_value.empty()) {
			return 0;
		} else {
			return 1;
		}
	}

	unsigned max_tokens() const {
		if (m_multitoken) {
			return 32000;
		} else if (m_zero_tokens) {
			return 0;
		} else {
			return 1;
		}
	}

	bool is_required() const { return m_required; }

	void xparse(boost::any& value_store, 
				const std::vector<std::basic_string<charT>>& new_tokens) const {
		if (new_tokens.empty() && !m_implicit_value.empty())
			value_store = m_implicit_value;
		else if (m_validate)
			m_validate(value_store, new_tokens);
		else
			validate(value_store, new_tokens, (T*)0, 0);
	}

	virtual bool apply_default(boost::any& value_store) const
	{
		if (m_default_value.empty()) {
			if (m_default_value_as_text.empty()) {
				return false;
			}
			else {
				std::vector<std::string> tokens;
				tokens.push_back(m_default_value_as_text);
				xparse(value_store, tokens);
				return true;
			}
		}
		else {
			value_store = m_default_value;
			return true;
		}
	}

	void notify(const boost::any& value_store) const
	{
		const T* value = boost::any_cast<T>(&value_store);
		if (m_store_to) {
			*m_store_to = *value;
		}
		if (m_notifier) {
			m_notifier(*value);
		}
	}

public: // typed_value_base overrides
	
	const std::type_info& value_type() const
	{
		return typeid(T);
	}
	

private:
	T* m_store_to;
	std::string m_value_name;
	boost::any m_default_value;
	std::string m_default_value_as_text;
	boost::any m_implicit_value;
	std::string m_implicit_value_as_text;
	bool m_composing, m_implicit, m_multitoken, m_zero_tokens, m_required;
	boost::function1<void, const T&> m_notifier;
	validate_function m_validate;
};

}}



template <class T>
boost::program_options::typed_value_extended<T>* value(T* v) {
	return new boost::program_options::typed_value_extended<T>(v);
}

template <class T>
boost::program_options::typed_value_extended<T>* value() {
	return value<T>(0);
}


#endif


