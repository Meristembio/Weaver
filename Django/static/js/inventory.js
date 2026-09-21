$(document).ready(function(){
    onReady();
});
function onReady(){
    $('.prevent_default').click(function(e){
        e.preventDefault();
    });
    /* dark mode */
    $('#btnSwitchDarkMode').unbind("click");
    $('#btnSwitchDarkMode').click(function(e){
        if (Cookies.get('dark_mode')  === 'true'){
            Cookies.set('dark_mode', false, { expires: 400 });
        }
        else{
            Cookies.set('dark_mode', true, { expires: 400 });
        }
        window.location.reload();
    });
    /* copy_clipboard */
    $('.copy_clipboard-child').click(function(e){
        e.preventDefault();
        copy_clipboard($(this).children('.copy_clipboard').first());
    });
    $('.copy_clipboard').click(function(e){
        e.preventDefault();
        copy_clipboard($(this));
    });
    function copy_clipboard(element){
        var $temp = $("<textarea>");
        $("body").append($temp);
        $temp.val(element.attr('data-cc')).select();
        document.execCommand("copy");
        $temp.remove();
        element.removeClass('bi-clipboard').addClass('bi-clipboard-check');
        setTimeout(function(){
            element.removeClass('bi-clipboard-check').addClass('bi-clipboard');
        }, 1000);
    }
    /* create/edit glycerolstock */
    if(typeof createglycerolstockPID !== 'undefined'){
    $('#id_plasmid').val(createglycerolstockPID);
    }
    $('.box-position-empty').click(function(){
        $('.box-position').removeClass('selected');
        $(this).addClass('selected');
        $('#id_box').val($(this).attr('data-box_id'));
        $('#id_box_column').val($(this).attr('data-box_column'));
        $('#id_box_row').val($(this).attr('data-box_row'));
    });
    /* select2 */
    var select2_ids = '#id_backbone, #id_inserts, #id_parent, #id_plasmid, #id_primer_f, #id_primer_r';
    $(select2_ids).select2();
    /* plasmid table massive actions */
    /* tooltips bootstrap */
    var tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'))
    var tooltipList = tooltipTriggerList.map(function (tooltipTriggerEl) {
      return new bootstrap.Tooltip(tooltipTriggerEl)
    });
    /* ove */
    if($("div#ove-viewer").length){
        $([document.documentElement, document.body]).animate({
            scrollTop: $("div#ove-viewer").offset().top
        }, 500);
    }
    /* digest */
    var digest_enzymes = [];
    $('#digest-choose button').click(function(){
        var enzymeName = $(this).attr('data-name');
        $(this).toggleClass('active');
        if(!digest_enzymes.includes(enzymeName)){
            digest_enzymes.push(enzymeName);
        } else{
            digest_enzymes.pop(enzymeName);
        }
        $('#digest_enzymes').val(JSON.stringify(digest_enzymes));
    });
    /* filters */
    var filter_buttons = $('.pe-table-filter-button');
    var filter_items = $('.filter-item');
    function resetTableFilters() {
        filter_buttons.removeClass('active');
        filter_items.removeClass('filter-hide');
    }
    filter_buttons.off("click").on("click", function(){
        filter_buttons.removeClass('active');
        $(this).addClass('active');
        var filter_id = $(this).attr('data-target');
        if(filter_id == 'all'){
            filter_items.removeClass('filter-hide');
            return;
        }
        filter_items.each(function(){
            var match = $(this).hasClass('filter-'+filter_id);
            if(!match){
                var find_in = $(this).children('td:first-child').children('a').attr('data-name');
                if(find_in)
                    match = find_in.startsWith(filter_id);
            }
            if(match){
                $(this).removeClass('filter-hide');
            } else{
                $(this).addClass('filter-hide');
            }
        });
    });
    $('#pe-table-filter-clear').off("click").on("click", function(){
        resetTableFilters();
        filter_buttons.filter('[data-target="all"]').first().addClass('active');
    });
    if (filter_buttons.length && !filter_buttons.filter('.active').length) {
        filter_buttons.filter('[data-target="all"]').first().addClass('active');
    }
    /* show_from_all_projects */
    $('#show_from_all_projects').click(function(){
        $(this).parent().submit();
    });
}

/* OVE save */
function saveOVE(sequenceDataToSave) {
    var genbankContents = window.bioParsers.jsonToGenbank(sequenceDataToSave, {isProtein: false});
    $.ajax({
      method: "POST",
      url: save_ove_path,
      data: {
        saveOve: true,
        gbContent: genbankContents
       },
      success: function(msg){
        window.toastr.success(msg.result);
      },
      error: function(XMLHttpRequest, textStatus, errorThrown) {
        window.toastr.error("Error while saving. Contact the administrator.");
      }
    })
}

function expandName(){
    var elements = $('#plasmids-table td a:first-child button.plasmid_list-name');
    if($('#table_search-expand').is(":checked")){
        elements.each(function(){
            $(this).attr('data-short', $(this).html())
            $(this).html($(this).parent().attr('data-name'))
        })
    } else {
        elements.each(function(){
            $(this).html($(this).attr('data-short'))
        })
    }
}
