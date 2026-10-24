document.addEventListener("DOMContentLoaded", function() {
    document.querySelectorAll('.highlight').forEach(block => {
        block.innerHTML = block.innerHTML.replace(/\bPASSED\b/g, '<span class="keyword-passed">PASSED</span>');
    });
});
