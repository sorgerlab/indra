// CURATION FUNCTIONS

// Force activate the sub items of the table of contents after page load
$(document).ready(function(){
    $('a[href="#statements"]').addClass('active')
})

// Variables
var pubmed_fetch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
var latestSubmission = {
                        'curator': '',
                        'api_key': '',
                        'ddSelect': '',
                        'ev_hash': '',
                        'source_hash': '',
                        'submit_status': 0
                    };

// Check the API key input
function keySubmit (key_value) {
    var ensure_user = document.getElementById("ensure_user_on_api_key")
    // Default value still there or nothing entered
    if (key_value == "No key given." | !key_value) {
        ensure_user.textContent = "No key given.";
    // nothing entered
    } else {
        console.log("document.getElementById(\"api_key_input\").value: " + document.getElementById("api_key_input").value);
        ensure_user.textContent = "Key stored!";
    }
}

function submitButtonClick(clickEvent) {
    // CURATOR
    var curator = document.getElementById("curator_input").value;
    if (!curator) {
        alert("Please enter a curator ID");
        return;
    }

    // API KEY
    var api_key = document.getElementById("api_key_input").value;
    if (!api_key) {
        alert("Please enter an API key!");
        return;
    }

    // Get mouseclick target, then parent's parent
    pn = clickEvent.target.parentNode.parentNode;
    btn_row_tag = pn.closest('tr');
    s = pn.getElementsByClassName("dropdown")[0]
        .getElementsByTagName("select")[0]

    // DROPDOWN SELECTION
    var err_select = s.options[s.selectedIndex].value;
    if (!err_select) {
        alert('Please select an error type or "correct" for the statement in the dropdown menu');
        return;
    }

    // TEXT BOX CONTENT
    // Get "form" node, then the value of "input"
    var user_text = pn
                .getElementsByClassName("form")[0]
                .getElementsByTagName("input")[0]
                .value;

    // Refuse submission if 'other' is selected without providing a description
    if (!user_text & err_select == "other") {
        alert('Must describe error when using option "other..."!')
        return;
    }

    // GET REFERENCE TO STATUS BOX (EMPTY UNTIL STATUS RECEIVED)
    var statusBox = pn.getElementsByClassName("submission_status")[0].children[0];

    // Step back to the preceding tr tag
    var pmid_row = btn_row_tag.previousElementSibling;

    // PMID
    // Get pmid_linktext content
    var pmid_text = pmid_row
                .getElementsByClassName("pmid_link")[0]
                .textContent.trim()

    // Icon 
    var icon = pmid_row.getElementsByClassName("curation_toggle")[0]

    // HASHES: source_hash & stmt_hash
    // source_hash == ev['source_hash'] == pmid_row.id; "evidence level"
    source_hash = pmid_row.id;
    // stmt_hash == hash == stmt_info['hash'] == table ID; "(pa-) statement level"
    stmt_hash = pmid_row.parentElement.parentElement.id;

    // CURATION DICT
    // example: curation_dict = {'tag': 'Reading', 'text': '"3200 A" is picked up as an agent.', 'curator': 'Klas', 'ev_hash': ev_hash};
    cur_dict = {
                'tag': err_select, 
                'text': user_text, 
                'curator': curator,
                'ev_hash': source_hash
            };

    // console.log("User: " + curator);
    // console.log("source hash: " + source_hash)
    // console.log("stmt hash: " + stmt_hash)
    // console.log("Error selected: " + err_select);
    // console.log("User feedback: " + user_text);
    // console.log("PMID: " + pmid_text);
    // console.log("cur_dict");
    // console.log(cur_dict);

    // SPAM CONTROL: preventing multiple clicks of the same curation in a row
    // If the new submission matches latest submission AND the latest submission was
    // successfully submitted, ignore the new submission
    if (latestSubmission['curator'] == curator &
        latestSubmission['api_key'] == document.getElementById("api_key_input").value &
        latestSubmission['ddSelect'] == err_select &
        latestSubmission['source_hash'] == source_hash &
        latestSubmission['stmt_hash'] == stmt_hash &
        latestSubmission['submit_status'] == 200) {
        alert('Already submitted curation successfully!');
        return;
    } else {
        latestSubmission['curator'] = curator;
        latestSubmission['api_key'] = document.getElementById("api_key_input").value;
        latestSubmission['ddSelect'] = err_select;
        latestSubmission['source_hash'] = source_hash;
        latestSubmission['stmt_hash'] = stmt_hash;
    };

    testing = false; // Set to true to test the curation endpoint of the API
    ajx_response = submitCuration(cur_dict, stmt_hash, statusBox, icon, testing);
    console.log("ajax response from submission: ");
    console.log(ajx_response);
};

// Submit curation
function submitCuration(curation_dict, hash, statusBox, icon, test) {
    // API KEY FOR ACCESSING DATABASE
    var api_key = document.getElementById("api_key_input").value;
    curation_addr = "https://db.indra.bio/curation/submit/";

    _url = curation_addr + hash + "?api_key=" + api_key;

    if (test) {
        console.log("Submitting test curation...");
        _url += "&test";
    };

    console.log("api key: " + api_key)
    console.log("url: " + _url);

    response = $.ajax({
        url: _url,
        type: "POST",
        dataType: "json",
        contentType: "application/json",
        data: JSON.stringify(curation_dict),
        complete: function(xhr, statusText) {
            latestSubmission['submit_status'] = xhr.status;
            switch (xhr.status) {
                case 200:
                    statusBox.textContent = "Curation submitted successfully!";
                    icon.style = "color: #00FF00" // Brightest green
                    break;
                case 400:
                    statusBox.textContent = xhr.status + ": Bad Curation Data";
                    icon.style = "color: #FF0000"; // Super red
                    break;
                case 404:
                    statusBox.textContent = xhr.status + ": Bad Link";
                    icon.style = "color: #FF0000";
                    break;
                case 500:
                    statusBox.textContent = xhr.status + ": Internal Server Error";
                    icon.style = "color: #FF0000";
                    break;
                case 504:
                    statusBox.textContent = xhr.status + ": Server Timeout";
                    icon.style = "color: #58D3F7"; // Icy blue
                    break;
                default:
                    console.log("Uncaught submission error: check ajax response")
                    console.log("xhr: ")
                    console.log(xhr)
                    statusBox.textContent = "Uncaught submission error; Code " + xhr.status;
                    icon.style = "color: #FF8000"; // Warning orange
                    break;
            }
        }
    });
    return response;
};

// Creates the dropdown div with the following structure
// <div class="dropdown" 
//      style="display:inline-block; vertical-align: middle;">
//     <select>
//         <option value="" selected disabled hidden>
//             Select error type...
//         </option>
//         <option value="correct">Correct</option>
//         <option value="entity_boundaries">Entity Boundaries</option>
//         <option value="grounding">Grounding</option>
//         <option value="no_relation">No Relation</option>
//         <option value="wrong_relation">Wrong Relation</option>
//         <option value="act_vs_amt">Activity vs. Amount</option>
//         <option value="polarity">Polarity</option>
//         <option value="negative_result">Negative Result</option>
//         <option value="hypothesis">Hypothesis</option>
//         <option value="agent_conditions">Agent Conditions</option>
//         <option value="mod_site">Modification Site</option>
//         <option value="other">Other...</option>
//     </select>
// </div>
function createDDDiv () {
    var ddContainer = document.createElement("div")
    ddContainer.className = "dropdown"
    ddContainer.style = "display:inline-block; vertical-align: middle; margin-left: 9%;"

    let ddSelect = document.createElement("select");

    // DROPDOWN OPTIONS
    // Default; This is the option placeholder
    placeholderOption = document.createElement("option");
    placeholderOption.value = "";
    placeholderOption.selected = "selected";
    placeholderOption.disabled = "disabled";
    placeholderOption.hidden = "hidden";
    placeholderOption.textContent = "Select error type..."
    ddSelect.appendChild(placeholderOption);
    // 1 "correct" No Error;
    option1 = document.createElement("option");
    option1.value = "correct";
    option1.textContent = "Correct";
    ddSelect.appendChild(option1);
    // 2 "entity_boundaries" Entity Boundaries;
    option2 = document.createElement("option");
    option2.value = "entity_boundaries";
    option2.textContent = "Entity Boundaries";
    ddSelect.appendChild(option2);
    // 3 "grounding" Grounding;
    option3 = document.createElement("option");
    option3.value = "grounding";
    option3.textContent = "Grounding";
    ddSelect.appendChild(option3);
    // 4 "no_relation" No Relation;
    option4 = document.createElement("option");
    option4.value = "no_relation";
    option4.textContent = "No Relation";
    ddSelect.appendChild(option4);
    // 5 "wrong_relation" Wrong Relation Type;
    option5 = document.createElement("option");
    option5.value = "wrong_relation";
    option5.textContent = "Wrong Relation";
    ddSelect.appendChild(option5);
    // 6 "act_vs_amt" Activity vs. Amount
    option6 = document.createElement("option");
    option6.value = "act_vs_amt";
    option6.textContent = "Activity vs. Amount";
    ddSelect.appendChild(option6);
    // 7 "polarity" Polarity;
    option7 = document.createElement("option");
    option7.value = "polarity";
    option7.textContent = "Polarity";
    ddSelect.appendChild(option7);
    // 8 "negative_result" Negative Result;
    option8 = document.createElement("option");
    option8.value = "negative_result";
    option8.textContent = "Negative Result";
    ddSelect.appendChild(option8);
    // 9 "hypothesis" Hypothesis;
    option9 = document.createElement("option");
    option9.value = "hypothesis";
    option9.textContent = "Hypothesis";
    ddSelect.appendChild(option9);
    // 10 "agent_conditions" Agent Conditions;
    option10 = document.createElement("option");
    option10.value = "agent_conditions";
    option10.textContent = "Agent Conditions";
    ddSelect.appendChild(option10);
    // 11 "mod_site" Modification Site;
    option11 = document.createElement("option");
    option11.value = "mod_site";
    option11.textContent = "Modification Site";
    ddSelect.appendChild(option11);
    // 12 "other" Other...
    option12 = document.createElement("option");
    option12.value = "other";
    option12.textContent = "Other...";
    ddSelect.appendChild(option12);
    // Add more options by following the structure above
    ddContainer.appendChild(ddSelect);
    return ddContainer;
};

// Creates the text box div with the following structure:
// <div class="form" 
//      style="display:inline-block; 
//             vertical-align: middle; 
//             top: 0px">
//     <form name="user_feedback_form">
//         <input type="text" 
//                maxlength="240"
//                name="user_feedback" 
//                placeholder="Optional description (240 chars)" 
//                value=""
//                style="width: 360px;">
//     </form>
// </div>
function createTBDiv () {
    var tbContainer = document.createElement("div")
    tbContainer.className = "form"
    tbContainer.style = "display:inline-block; vertical-align: middle; margin-left: 4%;"

    let tbForm = document.createElement("form")
    tbForm.name = "user_feedback_form"

    let tbInput = document.createElement("input")
    tbInput.type = "text";
    tbInput.maxlength = "240"
    tbInput.name = "user_feedback";
    tbInput.placeholder = "Optional description (240 chars)";
    tbInput.value = "";
    tbInput.style = "width: 360px;"

    tbForm.appendChild(tbInput);

    tbContainer.appendChild(tbForm);

    return tbContainer;
};

// Creates the submit button div with the following structure
// <div class="curation_button"
//      style="display:inline-block; 
//             vertical-align: middle;">
//     <button
//         type="button"
//         class="btn btn-default btn-submit pull-right"
//         style="padding: 2px 6px">Submit
//     </button>
//     < script tag type="text/javascript">
//     $(".btn-submit").off("click").on("click", function(b){
//         // Get parent node
//         parent_node = b.target.parentNode.parentNode
//         // Get reference to closest row tag (jquery)
//         this_row = $(this).closest("tr")
//         submitButtonClick(clickEvent)
//     })
//     </ script tag>
// </div>
function createSBDiv () {
    var sbContainer = document.createElement("div");
    sbContainer.className = "curation_button";
    sbContainer.style = "display:inline-block; vertical-align: middle; margin-left: 4%;";

    let sbButton = document.createElement("button");
    sbButton.type = "button";
    sbButton.className = "btn btn-default btn-submit pull-right";
    sbButton.style = "padding: 2px 6px; border: solid 1px #878787;";
    sbButton.textContent = "Submit";
    sbButton.onclick = submitButtonClick; // ATTACHES SCRIPT TO BUTTON

    sbContainer.appendChild(sbButton);

    return sbContainer;
};

// Creates the textbox that tells the user the status of the submission
// <div class="submission_status"
//      style="display:inline-block; 
//             vertical-align: middle;">
// <a class="submission_status"></a>
// </div>
function createStatusDiv () {
    var statusContainer = document.createElement("div");
    statusContainer.className = "submission_status";
    statusContainer.style = "display:inline-block; vertical-align: middle; margin-left: 4%;";

    let textContainer = document.createElement("i");
    textContainer.textContent = "";

    statusContainer.appendChild(textContainer);

    return statusContainer;

}

// Append row to the row that executed the click
// <tr class="cchild" style="border-top: 1px solid #FFFFFF;">
// <td colspan="4" style="padding: 0px; border-top: 1px solid #FFFFFF;">
function curationRowGenerator() {
    // Create new row element
    var newRow = document.createElement('tr');
    newRow.innerHTML = null;
    newRow.className = "cchild";
    newRow.style = "border-top: 1px solid #FFFFFF;";

    // Create new td element
    let newTD = document.createElement('td');
    newTD.style="padding: 0px; border-top: 1px solid #FFFFFF; white-space: nowrap; text-align: left;";
    newTD.setAttribute("colspan", "4");

    // Add dropdown div
    var dropdownDiv = createDDDiv();
    newTD.appendChild(dropdownDiv);
    // Add textbox 
    var textBoxDiv = createTBDiv();
    newTD.appendChild(textBoxDiv);
    // Add submit button
    var buttonDiv = createSBDiv();
    newTD.appendChild(buttonDiv);
    // Add submission response textbox
    var statusDiv = createStatusDiv();
    newTD.appendChild(statusDiv);

    // Add td to table row
    newRow.appendChild(newTD);

    return newRow;
};

// Adds the curation row to current
function addCurationRow(clickedRow) {
    // Generate new row
    var curationRow = curationRowGenerator();

    // Append new row to provided row
    clickedRow.parentNode.insertBefore(curationRow, clickedRow.nextSibling);
};

function getPubMedMETAxmlByPMID(pmid) {
    params_dict = {'db': 'pubmed',
        'retmode': 'xml',
        'rettype': 'docsum',
        'id': pmid
    };
    PubMedMETAxml = $.ajax({
        url: pubmed_fetch,
        type: "POST",
        dataType: "xml",
        data: params_dict,
    });
    return PubMedMETAxml
};

function pmidXML2dict(XML) {
    xml_dict = {};
    for (child of XML.children) {
        name = child.getAttribute("Name");
        type = child.getAttribute("Type");
        if (child.hasChildNodes() & type == "List") {
            // Javascript can't really do nice recursive functions...
            // special cases for "History" and "ArticleIds" which has unique inner Names
            if (name == "ArticleIds" | name == "History") {
                innerDict = {};
                for (c of child.children) {
                    innerDict[c.getAttribute("Name")] = c.textContent;
                }
                innerItems = innerDict;
            } else {
                innerList = [];
                for (c of child.children) {
                    innerList.push(c.textContent);
                }
                innerItems = innerList;
            }
            xml_dict[name] = innerItems
        } else if (child.tagName == "Item") {
            // Here just get the inner strings
            xml_dict[name] = child.textContent;
        } else if (child.tagName == "Id") {
            // Special case
            xml_dict["Id"] = child.textContent;
        } else {
            if (!xml_dict["no_key"]) {
                xml_dict["no_key"] = [child.textContent]
            } else {
                xml_dict["no_key"].push(child.textContent)
            }
        }
    }
    return xml_dict;
}

// Modify link hover text
function setPMIDlinkTitle(pmid, link_tag) {
    let pubmed_xml_promise = getPubMedMETAxmlByPMID(pmid);
    pubmed_xml_promise.then(function(responseXML) {
        docsum_xml = responseXML.getElementsByTagName('DocSum')[0];
        pmid_meta_dict = pmidXML2dict(docsum_xml);
        authorlist = pmid_meta_dict.AuthorList
        if (authorlist.length > 3) {
            authors = authorlist[0] + ", ... " + authorlist[authorlist.length-1];
        } else {
            authors = authorlist.join(", ");
        }
        // Shortened journal name is in .Source, while full name is in .FullJournalName
        journal = pmid_meta_dict.Source
        SO = pmid_meta_dict.SO
        title = pmid_meta_dict.Title

        pmid_hover_string = authors + ", \"" + title + "\", " + journal + ", " + SO
        link_tag.title = pmid_hover_string;
    })
}

// Loop all pmid link nodes and set title
function populatePMIDlinkTitles() {
    var pmid_link_array = document.getElementsByClassName("pmid_link");
    for (link_obj of pmid_link_array) {
        pmid = link_obj.textContent;
        setPMIDlinkTitle(pmid, link_obj)
    }
}

// Expand/collapse row
$(function() {
    $("td[class='curation_toggle']").click(function(event) {
        event.stopPropagation();
        var $target = $(event.target);
        if (event.target.dataset.clicked == "true") {
            // Toggle (animation duration in msec)
            $target.closest("tr").next().find("div").slideToggle(200);
        // First click event
        } else {
            // Stay down (animation duration in msec)
            $target.closest("tr").next().find("div").slideDown(400);

            // Change color of icon to light gray
            event.target.style="color:#A4A4A4;"

            // Set clicked to true
            event.target.dataset.clicked = "true"
        }
    });
});
